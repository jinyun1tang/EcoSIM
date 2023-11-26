module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use fileUtil, only : open_safe, check_read,int2str
  use minimathmod, only : isLeap
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use CanopyDataType
  use PlantTraitDataType
  use PlantMngmtDataType
  use EcosimConst
  use RootDataType
  use EcoSIMHistMod
  use CanopyRadDataType
  use GridDataType
  use EcoSIMCtrlMod
  use EcoSIMConfig

implicit none

  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  integer :: IDX,IMO,IYR,IDY,ICUT,IDYE,IDYG,IDYS,JCUT,LPY

  public :: ReadPlantInfo
  contains

  subroutine ReadPlantInfo(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!
! DESCRIPTION
!

  implicit none
  integer, intent(in) :: yearc, yeari
  integer, intent(in) :: NE   !
  integer, intent(in) :: NEX  !
  integer, intent(in) :: NHW,NHE,NVN,NVS     !simulation domain specification

! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
  if(lverb)WRITE(*,333)'ROUTQ'
  CALL ROUTQ(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  if(lverb)WRITE(*,333)'READQ'
  CALL READQ(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)

333   FORMAT(A8)

  end subroutine ReadPlantInfo
!------------------------------------------------------------------------------------------

  SUBROUTINE readq(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!!
! Description
! THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
!     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
!
  implicit none
  integer, intent(in) :: yearc,yeari,NEX
  integer, intent(in) :: NE,NHW,NHE,NVN,NVS
  integer :: NX,NY,NZ

! begin_execution

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      D9985: DO NZ=1,NP(NY,NX)
!
!       OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
!       FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!       each column has its own management
!       PREFIX=path for files in current or higher level directory
!       DATAP=PFT file name
        call ReadPlantProperties(NZ,NY,NX)

      ENDDO D9985
    ENDDO D9990
  ENDDO D9995

  call ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS)
  RETURN
  END SUBROUTINE readq
!------------------------------------------------------------------------------------------

  subroutine ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS)

  USE EcoSIMCtrlMod, only : pft_mgmt_in
  use netcdf
  use ncdio_pio
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS,yeari,yearc
  integer :: NY,NX,NZ

  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear,nyears,year
  integer :: NTOPO,ntopou
  integer :: M,NN,N,nn1
  real(r8) :: DY,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
  real(r8) :: ECUT24,HCUT,PCUT
  integer :: pft_nmgnt(JP)
  integer :: NH1,NV1,NH2,NV2,NS
  character(len=128) :: pft_pltinfo(JP),tstr
  character(len=128) :: pft_mgmtinfo(24,JP)

  if (len_trim(pft_mgmt_in)==0)then
    return
  else
  !  print*,'ReadPlantManagementNC'
  ! initialize the disturbance arrays
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO NZ=1,NP(NY,NX)
          DO M=1,366
            IHVST(NZ,M,NY,NX)=-1
            JHVST(NZ,M,NY,NX)=0
            HVST(NZ,M,NY,NX)=1.0E+06_r8
            THIN_pft(NZ,M,NY,NX)=-1.0_r8
            EHVST(1,ipld_leaf,NZ,M,NY,NX)=1.0_r8
            EHVST(1,ipld_nofoliar,NZ,M,NY,NX)=1.0_r8
            EHVST(1,ipld_woody,NZ,M,NY,NX)=1.0_r8
            EHVST(1,ipld_stdead,NZ,M,NY,NX)=1.0_r8
            EHVST(2,ipld_leaf,NZ,M,NY,NX)=1.0_r8
            EHVST(2,ipld_nofoliar,NZ,M,NY,NX)=1.0_r8
            EHVST(2,ipld_woody,NZ,M,NY,NX)=1.0_r8
            EHVST(2,ipld_stdead,NZ,M,NY,NX)=1.0_r8
          ENDDO
          PlantinDepth(NZ,NY,NX)=ppmc
        ENDDO
      ENDDO
    ENDDO

    call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

    call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))

    if(pft_dflag==0)then
      ! constant pft data
      iyear=1
      if(IGO>0)return
    else
      iyear=1
      DO while(.true.)
        call ncd_getvar(pftinfo_nfid,'year',iyear,year)
        if(year==yeari)exit
        iyear=iyear+1
      ENDDO
    endif
    !obtain the number of topographic units
    ntopou=get_dim_len(pftinfo_nfid, 'ntopou')
    nyears=get_dim_len(pftinfo_nfid, 'year')

    DO NTOPO=1,ntopou
      call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
      call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
      call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
      call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
      call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)

      call check_var(pftinfo_nfid, 'pft_pltinfo', vardesc, readvar)
      if(.not. readvar)then
        call endrun('fail to find pft_pltinfo in '//trim(mod_filename), __LINE__)
      endif

      call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_pltinfo, &
        start = (/1,1,NTOPO,iyear/),count = (/len(pft_pltinfo(1)),JP,1,1/)),&
        trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

      call check_var(pftinfo_nfid, 'nmgnts', vardesc, readvar)
      if(.not. readvar)then
        call endrun('fail to find pft_type in '//trim(mod_filename), __LINE__)
      endif

      call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_nmgnt, &
        start = (/1,NTOPO,iyear/),count = (/JP,1,1/)),trim(mod_filename)// &
        '::at line '//trim(int2str(__LINE__)))

      if(any(pft_nmgnt(1:NS)>0))then
        call check_var(pftinfo_nfid, 'pft_mgmt', vardesc, readvar)
        if(.not. readvar)then
          call endrun('fail to find pft_mgmt in '//trim(mod_filename), __LINE__)
        endif

        call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_mgmtinfo, &
          start = (/1,1,1,NTOPO,iyear/),count = (/len(pft_mgmtinfo(1,1)),24,JP,1,1/)),&
          trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
      endif

      DO NX=NH1,NH2
        DO NY=NV1,NV2
          DO NZ=1,MIN(NS,NP(NY,NX))
            tstr=trim(pft_pltinfo(NZ))
            read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
            read(tstr,*)DY,PPI(NZ,NY,NX),PlantinDepth(NZ,NY,NX)

            LPY=0
            if(isLeap(iyr) .and. IMO.GT.2)LPY=1
            !obtain the ordinal day
            IF(IMO.EQ.1)then
              IDY=IDX
            else
              IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
            endif
            IF(IDY.GT.0.AND.IYR.GT.0)THEN
              iDayPlanting(NZ,NY,NX)=IDY
              IYR=yearc
              iYearPlanting(NZ,NY,NX)=MIN(IYR,iYearCurrent)
              IDAYX(NZ,NY,NX)=iDayPlanting(NZ,NY,NX) !planting day
              IYRX(NZ,NY,NX)=iYearPlanting(NZ,NY,NX)   !planting year
              PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)     !population density
            ENDIF

            if(pft_nmgnt(NZ)>0)then
              NN=0
              DO nn1=1,pft_nmgnt(NZ)
                tstr=trim(pft_mgmtinfo(NN1,NZ))
                read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
                READ(TSTR,*)DY,ICUT,JCUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13,&
                   ECUT14,ECUT21,ECUT22,ECUT23,ECUT24

                if(isLeap(iyr) .and. IMO.GT.2)LPY=1
                !obtain the ordinal day
                IF(IMO.EQ.1)then
                  IDY=IDX
                else
                  IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
                endif

                IF(IDY.GT.0.AND.JCUT.EQ.1)THEN
                  iDayPlantHarvest(NZ,NY,NX)=IDY
                  IYR=yearc
                  iYearPlantHarvest(NZ,NY,NX)=MIN(IYR,iYearCurrent)
                ENDIF

                IHVST(NZ,IDY,NY,NX)=ICUT
                JHVST(NZ,IDY,NY,NX)=JCUT
                HVST(NZ,IDY,NY,NX)=HCUT
                THIN_pft(NZ,IDY,NY,NX)=PCUT
                EHVST(1,ipld_leaf,NZ,IDY,NY,NX)=ECUT11
                EHVST(1,ipld_nofoliar,NZ,IDY,NY,NX)=ECUT12
                EHVST(1,ipld_woody,NZ,IDY,NY,NX)=ECUT13
                EHVST(1,ipld_stdead,NZ,IDY,NY,NX)=ECUT14
                EHVST(2,ipld_leaf,NZ,IDY,NY,NX)=ECUT21
                EHVST(2,ipld_nofoliar,NZ,IDY,NY,NX)=ECUT22
                EHVST(2,ipld_woody,NZ,IDY,NY,NX)=ECUT23
                EHVST(2,ipld_stdead,NZ,IDY,NY,NX)=ECUT24

                IF(IHVST(NZ,IDY,NY,NX).EQ.4.OR.IHVST(NZ,IDY,NY,NX).EQ.6)THEN
                  !animal or insect biomass
                  NN=NN+1
                  if(mod(nn,2)==0)then
                    IDYE=IDY
                    D580: DO IDYG=IDYS+1,IDYE-1
                      IHVST(NZ,IDYG,NY,NX)=ICUT
                      JHVST(NZ,IDYG,NY,NX)=JCUT
                      HVST(NZ,IDYG,NY,NX)=HCUT
                      THIN_pft(NZ,IDYG,NY,NX)=PCUT
                      EHVST(1,ipld_leaf,NZ,IDYG,NY,NX)=ECUT11
                      EHVST(1,ipld_nofoliar,NZ,IDYG,NY,NX)=ECUT12
                      EHVST(1,ipld_woody,NZ,IDYG,NY,NX)=ECUT13
                      EHVST(1,ipld_stdead,NZ,IDYG,NY,NX)=ECUT14
                      EHVST(2,ipld_leaf,NZ,IDYG,NY,NX)=ECUT21
                      EHVST(2,ipld_nofoliar,NZ,IDYG,NY,NX)=ECUT22
                      EHVST(2,ipld_woody,NZ,IDYG,NY,NX)=ECUT23
                      EHVST(2,ipld_stdead,NZ,IDYG,NY,NX)=ECUT24
                    ENDDO D580
                  endif
                  IDYS=IDY
                ENDIF
              ENDDO
            endif
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    call ncd_pio_closefile(pftinfo_nfid)

  ENDIF

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP(NY,NX)
        IDAYY(NZ,NY,NX)=iDayPlantHarvest(NZ,NY,NX)
        IYRY(NZ,NY,NX)=iYearPlantHarvest(NZ,NY,NX)
      ENDDO
    ENDDO
  ENDDO
  end subroutine ReadPlantManagementNC
!------------------------------------------------------------------------------------------

  subroutine ReadPlantProperties(NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX

  integer :: N, NB
  real(r8) :: VRNXI,VRNLI

! begin_execution

!
! READ INPUTS FOR EACH PLANT SPECIES
!
  IF(DATAP(NZ,NY,NX).NE.'NO')THEN
    write(101,*)'NZ=',NZ,DATAP(NZ,NY,NX)
    call ReadPlantTraitsNC(NZ,NY,NX,VRNLI,VRNXI)
!
!   RE-CALCULATE PLANT INPUTS IN MODEL UNITS
!
!   CanopySWabsorpty_pft,CanopyPARabsorpty_pft=leaf SW,PAR absorbtivity
!
    VCMX(NZ,NY,NX)=2.5_r8*VCMX(NZ,NY,NX)
    VOMX(NZ,NY,NX)=2.5_r8*VOMX(NZ,NY,NX)
    VCMX4(NZ,NY,NX)=2.5_r8*VCMX4(NZ,NY,NX)
    ETMX(NZ,NY,NX)=2.5_r8*ETMX(NZ,NY,NX)
    CanopySWabsorpty_pft(NZ,NY,NX)=1.0_r8-CanopySWAlbedo_pft(NZ,NY,NX)-TAUR(NZ,NY,NX)
    CanopyPARabsorpty_pft(NZ,NY,NX)=1.0_r8-CanopyPARalbedo_pft(NZ,NY,NX)-TAUP(NZ,NY,NX)
    CanopySWAlbedo_pft(NZ,NY,NX)=CanopySWAlbedo_pft(NZ,NY,NX)/CanopySWabsorpty_pft(NZ,NY,NX)
    CanopyPARalbedo_pft(NZ,NY,NX)=CanopyPARalbedo_pft(NZ,NY,NX)/CanopyPARabsorpty_pft(NZ,NY,NX)
    TAUR(NZ,NY,NX)=TAUR(NZ,NY,NX)/CanopySWabsorpty_pft(NZ,NY,NX)
    TAUP(NZ,NY,NX)=TAUP(NZ,NY,NX)/CanopyPARabsorpty_pft(NZ,NY,NX)
    ANGBR(NZ,NY,NX)=SIN(ANGBR(NZ,NY,NX)/57.29578_r8)
    ANGSH(NZ,NY,NX)=SIN(ANGSH(NZ,NY,NX)/57.29578_r8)
    MatureGroup_pft(NZ,NY,NX)=GROUPX(NZ,NY,NX)

    IF(iPlantTurnoverPattern(NZ,NY,NX).NE.0)THEN
!
      RefNodeInitRate(NZ,NY,NX)=RefNodeInitRate(NZ,NY,NX)/MaxNodesPerBranch
      RefLeafAppearRate(NZ,NY,NX)=RefLeafAppearRate(NZ,NY,NX)/MaxNodesPerBranch
      MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)/MaxNodesPerBranch
      XTLI(NZ,NY,NX)=XTLI(NZ,NY,NX)/MaxNodesPerBranch
    ENDIF
    MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)-XTLI(NZ,NY,NX)
    IF(XDL(NZ,NY,NX).LT.0.0_r8)THEN
      XDL(NZ,NY,NX)=DayLenthMax(NY,NX)
    ENDIF
    D5: DO NB=1,JC
      IF(iPlantPhenologyType(NZ,NY,NX).EQ.0.AND.iPlantPhenologyPattern(NZ,NY,NX).NE.iplt_annual)THEN
        HourThreshold4LeafOut(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNLI+144.0_r8*iPlantInitThermoAdaptZone(NZ,NY,NX)*(NB-1))
        HourThreshold4LeafOff(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNXI+144.0_r8*iPlantInitThermoAdaptZone(NZ,NY,NX)*(NB-1))
      ELSE
        HourThreshold4LeafOut(NB,NZ,NY,NX)=VRNLI
        HourThreshold4LeafOff(NB,NZ,NY,NX)=VRNXI
      ENDIF
    ENDDO D5
  ENDIF
! WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,XDL(NZ,NY,NX)
!1111    FORMAT(A20,2I8,E12.4)
  end subroutine ReadPlantProperties

!------------------------------------------------------------------------------------------
  subroutine ReadPlantTraitsNC(NZ,NY,NX,VRNLI,VRNXI)
  use EcoSIMCtrlMod, only : pft_nfid
  use GrosubPars, only : get_pft_loc
  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(out) :: VRNLI,VRNXI
  integer :: loc,N

  loc=get_pft_loc(DATAP(NZ,NY,NX)(1:6))
  DATAPI(NZ,NY,NX)=loc
  call ncd_getvar(pft_nfid, 'ICTYP', loc, iPlantPhotosynthesisType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IGTYP', loc, iPlantMorphologyType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ISTYP', loc, iPlantPhenologyPattern(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IDTYP', loc, iPlantDevelopPattern(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'INTYP', loc, iPlantNfixType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IWTYP', loc, iPlantPhenologyType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IPTYP', loc, iPlantPhotoperiodType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IBTYP', loc, iPlantTurnoverPattern(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IRTYP', loc, iPlantGrainType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'MY', loc, MY(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ZTYPI', loc, iPlantInitThermoAdaptZone(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'VCMX', loc, VCMX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VOMX', loc, VOMX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VCMX4', loc, VCMX4(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO2', loc, XKCO2(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKO2', loc, XKO2(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO24', loc, XKCO24(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RUBP', loc, RUBP(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PEPC', loc, PEPC(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ETMX', loc, ETMX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL', loc, CHL(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL4', loc, CHL4(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'FCO2', loc, CanPCi2CaRatio(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'ALBR', loc, CanopySWAlbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ALBP', loc, CanopyPARalbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUR', loc, TAUR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUP', loc, TAUP(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'XRNI', loc, RefNodeInitRate(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XRLA', loc, RefLeafAppearRate(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CTC', loc, TCelciusChill4Seed(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VRNLI', loc,VRNLI)
  call ncd_getvar(pft_nfid, 'VRNXI', loc,VRNXI)
  call ncd_getvar(pft_nfid, 'WDLF', loc,WDLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PB', loc,MinNonstructalC4InitBranch(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'GROUPX', loc,GROUPX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XTLI', loc,XTLI(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XDL', loc,XDL(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XPPD', loc,XPPD(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'SLA1', loc,SLA1(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SSL1', loc,SSL1(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SNL1', loc,SNL1(NZ,NY,NX))


  call ncd_getvar(pft_nfid, 'CLASS', loc,CLASS(1:JLI,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CFI', loc,ClumpFactort0(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGBR', loc,ANGBR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGSH', loc,ANGSH(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'STMX', loc,STMX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SDMX', loc,SDMX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRMX', loc,MaxSeedCMass(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRDM', loc,SeedCMass(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GFILL', loc,GFILL(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'WTSTDI', loc,StandingDeadInitC_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'RRAD1M', loc,MaxPrimRootRadius(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RRAD2M', loc,MaxSecndRootRadius(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PORT', loc,RootPorosity(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PR', loc,MinNonstructuralC4InitRoot(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRR', loc,RSRR(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRA', loc,RSRA(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PTSHT', loc,ShutRutNonstructElmntConducts(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RTFQ', loc,RTFQ(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZH', loc,UPMXZH(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZH', loc,UPKMZH(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZH', loc,UPMNZH(1,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZO', loc,UPMXZO(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZO', loc,UPKMZO(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZO', loc,UPMNZO(1,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXPO', loc,UPMXPO(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMPO', loc,UPKMPO(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNPO', loc,UPMNPO(1,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'OSMO', loc,OSMO(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RCS', loc,RCS(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSMX', loc,RSMX(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'DMLF', loc,LeafBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSHE', loc,PetioleBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSTK', loc,StalkBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRSV', loc,ReserveBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMHSK', loc,HuskBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMEAR', loc,EarBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMGR', loc,GrainBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRT', loc,RootBiomGrowthYield(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMND', loc,DMND(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CNLF', loc,CNLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSHE', loc,CNSHE(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSTK', loc,rNCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRSV', loc,CNRSV(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNHSK', loc,CNHSK(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNEAR', loc,CNEAR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNGR', loc,CNGR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNumRootAxes_pft', loc,RootrNC_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNND', loc,CNND(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CPLF', loc,CPLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSHE', loc,CPSHE(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSTK', loc,CPSTK(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRSV', loc,CPRSV(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPHSK', loc,CPHSK(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPEAR', loc,CPEAR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPGR', loc,CPGR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRT', loc,CPRT(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPND', loc,CPND(NZ,NY,NX))

  if(lverb)then
    call pft_display(NZ,NY,NX)
    call photosyns_trait_disp(NZ,NY,NX)
    call plant_optic_trait_disp(NZ,NY,NX)
    call Phenology_trait_disp(NZ,NY,NX,VRNLI,VRNXI)
    call morphology_trait_disp(NZ,NY,NX)
    call Root_trait_disp(NZ,NY,NX)
    call Root_nutrient_trait_disp(NZ,NY,NX)
    call plant_water_trait_disp(NZ,NY,NX)
    call plant_biomyield_trait_disp(NZ,NY,NX)
    call plant_biomstoich_trait_disp(NZ,NY,NX)
  endif

  end subroutine ReadPlantTraitsNC
!------------------------------------------------------------------------------------------
  subroutine pft_display(NZ,NY,NX)
  use abortutils , only : endrun
  implicit none
  integer, intent(in) :: NZ,NY,NX
!   iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4
!   iPlantMorphologyType=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   iPlantPhenologyPattern=growth habit:0=annual,1=perennial
!   iPlantDevelopPattern=growth habit:0=determinate,1=indetermimate
!   iPlantNfixType=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   iPlantPhenologyType=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   iPlantPhotoperiodType=photoperiod type:0=day neutral,1=short day,2=long day
!   iPlantTurnoverPattern=turnover:if iPlantMorphologyType=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if iPlantMorphologyType=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   iPlantGrainType=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   iPlantInitThermoAdaptZone=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical

  write(*,*)'PLANT FUNCTIONAL TYPE for (NZ,NY,NX)=',NZ,NY,NX
  select CASE (iPlantPhotosynthesisType(NZ,NY,NX))
  case (3)
    write(*,*)'C3 photosynthesis'
  case (4)
    write(*,*)'C4 photosynthesis'
  case default
    write(*,*)'photosynthesis type not defined'
  end select

  select case(iPlantMorphologyType(NZ,NY,NX))
  case (0)
    write(*,*)'shallow root profile, like bryophytes'
  case (1)
    write(*,*)'intermediate root profile, like herbs'
  case (2)
    write(*,*)'deep root profile, like trees'
  case default
    write(*,*)'root profile not defined'
  end select

  select case (iPlantPhenologyPattern(NZ,NY,NX))
  case (0)
    write(*,*)'Annual plant'
  case (1)
    write(*,*)'perennial plant'
  case default
    write(*,*)'growth habit not defined'
  end select

  select case (iPlantDevelopPattern(NZ,NY,NX))
  case (0)
    write(*,*)'determinate growth habit'
  case (1)
    write(*,*)'indetermimate growth habit'
  case default
    write(*,*)'growth habit not defined'
  end select

  select case (iPlantNfixType(NZ,NY,NX))
! 1,2, 3, e.g. legumes
  case (1)
    write(*,*)'Rapid root N2 fixation symbiosis'
  case (2)
    write(*,*)'Intermediate root N2 fixation symbiosis'
  case (3)
    write(*,*)'Slow root N2 fixation symbiosis'
!4,5,6, e.g. cyanobacteria
  case (4)
    write(*,*)'Rapid canopy N2 fixation symbiosis'
  case (5)
    write(*,*)'Intermediate canopy N2 fixation symbiosis'
  case (6)
    write(*,*)'Slow canopy N2 fixation symbiosis'
  case default
    write(*,*)'No N2 fixation symbiosis defined'
  end select

  select case(iPlantPhenologyType(NZ,NY,NX))
  case (0)
    write(*,*)'phenology type: evergreen'
  case (1)
    write(*,*)'phenology type: cold deciduous'
  case (2)
    write(*,*)'phenology type: drought deciduous'
  case (3)
    write(*,*)'phenology type: cold+drought deciduous'
  case default
    write(*,*)'phenology type not defined'
  end select

  select case(iPlantPhotoperiodType(NZ,NY,NX))
  case (0)
    write(*,*)'day neutral photoperiod'
  case (1)
    write(*,*)'short day photoperiod'
  case (2)
    write(*,*)'long day photoperiod'
  case default
    write(*,*)'photoperiod not defined'
  end select

  select case(iPlantTurnoverPattern(NZ,NY,NX))
  case (0, 1)
    write(*,*)'plant biome turnover rapid (deciduous)'
  case (2)
    write(*,*)'Plant biome turnover very slow (needleleaf evergreen)'
  case (3)
    write(*,*)'Plant biome turnover slow broadleaf evergreen'
  case (4)
    write(*,*)'Plant biome turnover semi-deciduous'
  case (5)
    write(*,*)'Plant biome turnover semi-evergreen'
  case default
    write(*,*)'Plant biome turnover not defined'
  end select

  select case(iPlantGrainType(NZ,NY,NX))
  case (0)
    write(*,*)'Above ground storage organ'
  case (1)
    write(*,*)'Belowground storage organ'
  case default
    write(*,*)'Storage organ not defined'
  end select

  select case(MY(NZ,NY,NX))
  case (1)
    write(*,*)'No mycorrhizal'
  case (2)
    write(*,*)'Mycorrhizal'
  case default
    write(*,*)'Wrong option for mycorrhizae'
  end select

  select case(INT(iPlantInitThermoAdaptZone(NZ,NY,NX)))
  case (1)
    write(*,*)'thermal adaptation zone: arctic, boreal'
  case (2)
    write(*,*)'thermal adaptation zone: cool temperate'
  case (3)
    write(*,*)'thermal adaptation zone: warm temperate'
  case (4)
    write(*,*)'thermal adaptation zone: subtropical'
  case (5)
    write(*,*)'thermal adaptation zone: tropical'
  case default
    write(*,*)'Not thermal adaptation zone defined'
  end select

  end subroutine pft_display

!------------------------------------------------------------------------------------------
  subroutine photosyns_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'PHOTOSYNTHETIC PROPERTIES'
  write(*,*)'Specific rubisco carboxylase (umol C g-1 s-1): ',VCMX(NZ,NY,NX)
  write(*,*)'Specific rubisco oxygenase (umol O g-1 s-1): ',VOMX(NZ,NY,NX)
  write(*,*)'Specific PEP carboxylase activity (umol g-1 s-1): '&
    ,VCMX4(NZ,NY,NX)
  write(*,*)'Km for VCMX (uM): ',XKCO2(NZ,NY,NX)
  write(*,*)'Km for VOMX (uM): ',XKO2(NZ,NY,NX)
  write(*,*)'KM for VCMX4 (uM): ',XKCO24(NZ,NY,NX)
  write(*,*)'Fraction of leaf protein in rubisco: ',RUBP(NZ,NY,NX)
  write(*,*)'Fraction of leaf protein in PEP carboxylase: ',PEPC(NZ,NY,NX)
  write(*,*)'Specific chlorophyll activity (umol e- g-1 s-1): ',ETMX(NZ,NY,NX)
  if(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic3_photo)then
    write(*,*)'Fraction of leaf protein as chlorophyll in mesophyll (C3) ',CHL(NZ,NY,NX)
  elseif(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic4_photo)then
    write(*,*)'Fraction of leaf protein as chlorophyll in bundle sheath(C4)',CHL(NZ,NY,NX)
    write(*,*)'fraction of leaf protein in mesophyll chlorophyll(C4)',CHL4(NZ,NY,NX)
  endif
  write(*,*)'intercellular:atmospheric CO2 concentration ratio:',CanPCi2CaRatio(NZ,NY,NX)
  end subroutine photosyns_trait_disp

!------------------------------------------------------------------------------------------

  subroutine Phenology_trait_disp(NZ,NY,NX,VRNLI,VRNXI)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: VRNLI,VRNXI
  write(*,*)'PHENOLOGICAL PROPERTIES'
  write(*,*)'rate of node initiation at 25oC (h-1): XRNI',RefNodeInitRate(NZ,NY,NX)
  write(*,*)'rate of leaf appearance at 25oC (h-1): XRLA',RefLeafAppearRate(NZ,NY,NX)
  write(*,*)'chilling temperature for CO2 fixation, '// &
    'seed loss (oC): CTC',TCelciusChill4Seed(NZ,NY,NX)
  write(*,*)'hour requirement for spring leafout: VRNLI',VRNLI
  write(*,*)'hour requirement for autumn leafoff: VRNXI',VRNXI
  write(*,*)'leaf length:width ratio: WDLF',WDLF(NZ,NY,NX)
  write(*,*)'nonstructural C concentration needed for branching'// &
    ':PB',MinNonstructalC4InitBranch(NZ,NY,NX)
  end subroutine Phenology_trait_disp
!------------------------------------------------------------------------------------------
  subroutine morphology_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer :: N

  write(*,*)'MORPHOLOGICAL PROPERTIES'
  write(*,*)'growth in leaf area vs mass: SLA1 ',SLA1(NZ,NY,NX)
  write(*,*)'growth in petiole length vs mass: SSL1 ',SSL1(NZ,NY,NX)
  write(*,*)'growth in internode length vs mass: SNL1',SNL1(NZ,NY,NX)
  write(*,*)'fraction of leaf area in 0-22.5,45,67.5,90o '// &
    'inclination classes: CLASS',(CLASS(N,NZ,NY,NX),N=1,JLI)
  write(*,*)'initial clumping factor: CFI',ClumpFactort0(NZ,NY,NX)
  write(*,*)'stem angle from horizontal: ANGBR',ANGBR(NZ,NY,NX)
  write(*,*)'petiole angle from horizontal: ANGSH',ANGSH(NZ,NY,NX)
  write(*,*)'maximum potential seed mumber from '// &
    'pre-anthesis stalk growth: STMX',STMX(NZ,NY,NX)
  write(*,*)'maximum seed number per STMX: SDMX',SDMX(NZ,NY,NX)
  write(*,*)'maximum seed size per SDMX (g): GRMX',MaxSeedCMass(NZ,NY,NX)
  write(*,*)'seed size at planting (g): GRDM',SeedCMass(NZ,NY,NX)    !could be greater than MaxSeedCMass, accouting for seedling
  write(*,*)'grain filling rate at 25 oC (g seed-1 h-1): GFILL',GFILL(NZ,NY,NX)
  write(*,*)'mass of dead standing biomass at planting: WTSTDI',StandingDeadInitC_pft(NZ,NY,NX)
  end subroutine morphology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'ROOT CHARACTERISTICS'
  write(*,*)'radius of primary roots: MaxPrimRootRadius',MaxPrimRootRadius(1,NZ,NY,NX)
  write(*,*)'radius of secondary roots: MaxSecndRootRadius',MaxSecndRootRadius(1,NZ,NY,NX)
  write(*,*)'primary/fine root porosity: PORT',RootPorosity(1,NZ,NY,NX)
  write(*,*)'nonstructural C concentration needed for root'// &
    ' branching: PR',MinNonstructuralC4InitRoot(NZ,NY,NX)
  write(*,*)'radial root resistivity for water uptake (m2 MPa-1 h-1): RSRR',RSRR(1,NZ,NY,NX)
  write(*,*)'axial root resistivity for water uptake (m2 MPa-1 h-1): RSRA',RSRA(1,NZ,NY,NX)
  write(*,*)'rate constant for equilibrating shoot-root '// &
    'nonstructural C concn: PTSHT',ShutRutNonstructElmntConducts(NZ,NY,NX)
  write(*,*)'root branching frequency (m-1): RTFQ',RTFQ(NZ,NY,NX)
  end subroutine Root_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_nutrient_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'ROOT UPTAKE PARAMETERS'
  write(*,*)'NH4 max uptake (g m-2 h-1): UPMXZH',UPMXZH(1,NZ,NY,NX)
  write(*,*)'NH4 uptake Km (uM): UPKMZH',UPKMZH(1,NZ,NY,NX)
  write(*,*)'NH4 uptake min concn (uM): UPMNZH',UPMNZH(1,NZ,NY,NX)
  write(*,*)'NO3 max uptake (g m-2 h-1): UPMXZO',UPMXZO(1,NZ,NY,NX)
  write(*,*)'NO3 uptake Km (uM): UPKMZO',UPKMZO(1,NZ,NY,NX)
  write(*,*)'NO3 uptake min concn (uM): UPMNZO',UPMNZO(1,NZ,NY,NX)
  write(*,*)'H2PO4 max uptake (g m-2 h-1): UPMXPO',UPMXPO(1,NZ,NY,NX)
  write(*,*)'H2PO4 uptake Km (uM): UPKMPO',UPKMPO(1,NZ,NY,NX)
  write(*,*)'H2PO4 uptake min conc (uM): UPMNPO',UPMNPO(1,NZ,NY,NX)
  end subroutine Root_nutrient_trait_disp


!------------------------------------------------------------------------------------------
  subroutine plant_water_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'WATER RELATIONS'
  write(*,*)'leaf osmotic potential at zero leaf water '// &
     'potential (MPa): OSMO',OSMO(NZ,NY,NX)
  write(*,*)'shape parameter for stomatal resistance vs '// &
     'leaf turgor potential: RCS',RCS(NZ,NY,NX)
  write(*,*)'cuticular resistance (s m-1): RSMX',RSMX(NZ,NY,NX)
  end subroutine plant_water_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomyield_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX


  write(*,*)'ORGAN GROWTH YIELDS'
  write(*,*)'leaf dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMLF',LeafBiomGrowthYield(NZ,NY,NX)
  write(*,*)'petiole dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMSHE',PetioleBiomGrowthYield(NZ,NY,NX)
  write(*,*)'stalk dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMSTK',StalkBiomGrowthYield(NZ,NY,NX)
  write(*,*)'stalk reserve C production vs '// &
    'nonstructural C consumption (g g-1): DMRSV',ReserveBiomGrowthYield(NZ,NY,NX)
  write(*,*)'husk dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1): DMHSK',HuskBiomGrowthYield(NZ,NY,NX)
  write(*,*)'ear dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1): DMEAR',EarBiomGrowthYield(NZ,NY,NX)
  write(*,*)'grain C production vs nonstructural C'// &
    ' consumption (g g-1): DMGR',GrainBiomGrowthYield(NZ,NY,NX)
  write(*,*)'root dry matter C production vs nonstructural C'// &
    ' consumption (g g-1): DMRT',RootBiomGrowthYield(NZ,NY,NX)
  write(*,*)'nodule bacteria in root nodule,canopy dry matter'// &
    'C production vs nonstructural C consumption (g g-1): DMND' &
    ,DMND(NZ,NY,NX)
  end subroutine plant_biomyield_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomstoich_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'ORGAN N AND P CONCENTRATIONS'
  write(*,*)'NC ratio in plant leaves: CNLF',CNLF(NZ,NY,NX)
  write(*,*)'NC ratio in plant petiole: CNSHE',CNSHE(NZ,NY,NX)
  write(*,*)'NC ratio in plant stalk: CNSTK',rNCStalk_pft(NZ,NY,NX)
  write(*,*)'NC ratio in plant stalk reserve: CNRSV',CNRSV(NZ,NY,NX)
  write(*,*)'NC ratio in plant husk: CNHSK',CNHSK(NZ,NY,NX)
  write(*,*)'NC ratio in plant ear: CNEAR',CNEAR(NZ,NY,NX)
  write(*,*)'NC ratio in plant grain: CNGR',CNGR(NZ,NY,NX)
  write(*,*)'NC ratio in plant root: CNumRootAxes_pft',RootrNC_pft(NZ,NY,NX)
  write(*,*)'NC ratio in plant nodule: CNND',CNND(NZ,NY,NX)
  write(*,*)'PC ratio in plant leaves: CPLF',CPLF(NZ,NY,NX)
  write(*,*)'PC ratio in plant petiole: CPSHE',CPSHE(NZ,NY,NX)
  write(*,*)'PC ratio in plant stalk: CPSTK',CPSTK(NZ,NY,NX)
  write(*,*)'PC ratio in plant stalk reserve: CPRSV',CPRSV(NZ,NY,NX)
  write(*,*)'PC ratio in plant husk: CPHSK',CPHSK(NZ,NY,NX)
  write(*,*)'PC ratio in plant ear: CPEAR',CPEAR(NZ,NY,NX)
  write(*,*)'PC ratio in plant grain: CPGR',CPGR(NZ,NY,NX)
  write(*,*)'PC ratio in plant root: CPRT',CPRT(NZ,NY,NX)
  write(*,*)'PC ratio in plant nodule: CPND',CPND(NZ,NY,NX)
  end subroutine plant_biomstoich_trait_disp

!------------------------------------------------------------------------------------------

  subroutine plant_optic_trait_disp(NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX

  write(*,*)'OPTICAL PROPERTIES'
  write(*,*)'leaf SW albedo: ALBR',CanopySWAlbedo_pft(NZ,NY,NX)
  write(*,*)'leaf PAR albedo: ALBP',CanopyPARalbedo_pft(NZ,NY,NX)
  write(*,*)'leaf SW transmission: TAUR',TAUR(NZ,NY,NX)
  write(*,*)'leaf PAR transmission: TAUP',TAUP(NZ,NY,NX)
  end subroutine plant_optic_trait_disp
!------------------------------------------------------------------------------------------

  SUBROUTINE routq(yearc,yeari,NE,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE OPENS CHECKPOINT FILES AND READS
!     FILE NAMES FOR PLANT SPECIES AND MANAGEMENT
!
      use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  integer, intent(in) :: yearc,yeari
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS

  CHARACTER(len=16) :: OUTX,OUTC,OUTM,OUTR,OUTQ
  CHARACTER(len=4) :: CHARY
  integer :: IDATE,NY,NX,NZ

! begin_execution
!
! OPEN CHECKPOINT FILES FOR PLANT VARIABLES
!
  IF(is_first_year)THEN
    D9999: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        DO NZ=1,JP
          iDayPlanting(NZ,NY,NX)=-1E+06
          iYearPlanting(NZ,NY,NX)=-1E+06
          iDayPlantHarvest(NZ,NY,NX)=1E+06
          iYearPlantHarvest(NZ,NY,NX)=1E+06
        enddo
      enddo
    ENDDO D9999

  ENDIF
  if(lverb)write(*,*)'ReadPlantInfoNC'
  call ReadPlantInfoNC(yeari,NE,NEX,NHW,NHE,NVN,NVS)

  END subroutine routq


!------------------------------------------------------------------------------------------


  subroutine ReadPlantInfoNC(yeari,NE,NEX,NHW,NHE,NVN,NVS)

  USE EcoSIMCtrlMod, only : pft_mgmt_in
  use netcdf
  use ncdio_pio
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: NE,NEX,NHW,NHE,NVN,NVS
  integer, intent(in) :: yeari
  integer :: IDATE
  integer :: NPP(JY,JX)
  integer :: IYR,NX,NY,NZ,NN,NH1,NH2,NV1,NV2,NS
  CHARACTER(len=2) :: CLIMATE
  integer :: ntopou,NTOPO
  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear,year
  integer :: ret
  character(len=10) :: pft_gtype(JP)

  if (len_trim(pft_mgmt_in)==0)then
   !no plant management info to read
    D5995: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NP(NY,NX)=0
        DO NZ=1,NP0(NY,NX)
          DATAP(NZ,NY,NX)='NO'
          DATAM(NZ,NY,NX)='NO'
        enddo
      enddo
    ENDDO D5995
  else

    call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

    call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)

    if(.not. readvar)then
      call endrun('fail to find pft_dflag in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
      trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

    if (pft_dflag==-1)then
      !no pft data
      DO NX=NHW,NHE
        DO  NY=NVN,NVS
          NP(NY,NX)=0
          DO NZ=1,NP0(NY,NX)
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          enddo
        enddo
      ENDDO
    else
      if(pft_dflag==0)then
        ! constant pft data
        iyear=1
      else
        iyear=1
        DO while(.true.)
          call ncd_getvar(pftinfo_nfid,'year',iyear,year)
          if(year==yeari)exit
          iyear=iyear+1
        ENDDO
      endif
      !obtain the number of topographic units
      ntopou=get_dim_len(pftinfo_nfid, 'ntopou')

      DO NTOPO=1,ntopou
        call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
        call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
        call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
        call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)

        NN=1
        D4995: DO NX=NH1,NH2
          D4990: DO NY=NV1,NV2
            NP0(NY,NX)=NS
            D4985: DO NZ=1,NS
              LSG(NZ,NY,NX)=NN
            ENDDO D4985
          ENDDO D4990
        ENDDO D4995

        IF(NS.GT.0)THEN
          ! a vegetated topo unit
          call check_var(pftinfo_nfid, 'pft_type', vardesc, readvar)
          if(.not. readvar)then
            call endrun('fail to find pft_type in '//trim(mod_filename), __LINE__)
          endif
          call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_gtype, &
            start = (/1,1,NTOPO,iyear/),count = (/len(pft_gtype(1)),JP,1,1/)), &
            trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

          D4975: DO NX=NH1,NH2
            D4970: DO NY=NV1,NV2
              D4965: DO NZ=1,NS
                IF(IETYP(NY,NX).GT.0)THEN
                  WRITE(CLIMATE,'(I2)')IETYP(NY,NX)
                  !the type of pft is specified by genra+Koppen climate zone
                  DATAX(NZ)=pft_gtype(NZ)(1:4)//CLIMATE
                ENDIF

              ENDDO D4965
            ENDDO D4970
          ENDDO D4975
        ENDIF

!        IF(.not. is_restart())THEN
    ! there was no chechk point file read in, so update pft info
    ! from input file
          D8995: DO NX=NH1,NH2
            D8990: DO NY=NV1,NV2
              NP(NY,NX)=NS
    !DATAP(NZ,NY,NX) and DATAM(NZ,NY,NX) are to be read in readqmod.F90
              D100: DO NZ=1,NP(NY,NX)
                DATAP(NZ,NY,NX)=DATAX(NZ)
              ENDDO D100

              D101: DO NZ=NP(NY,NX)+1,JP
                DATAP(NZ,NY,NX)='NO'
              ENDDO D101
            ENDDO D8990
          ENDDO D8995
!        ELSE
!read from chck point file, i.e. datap and datam
!          call read_checkpt(NS,NH1,NH2,NV1,NV2,NHW,NHE,NVN,NVS)
!        ENDIF
      ENDDO
    ENDIF
    call ncd_pio_closefile(pftinfo_nfid)
  endif
  end subroutine ReadPlantInfoNC

!------------------------------------------------------------------------------------------

  subroutine read_checkpt(NS,NH1,NH2,NV1,NV2,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NS,NH1,NH2,NV1,NV2,NHW,NHE,NVN,NVS
  integer :: NY,NX,NZ,NN
  integer :: IDATE,IYR
  integer :: NPP(JY,JX)
  CHARACTER(len=16) :: DATAA(JP,JY,JX),DATAB(JP,JY,JX)

  NPP =0

  REWIND(30)

8000  CONTINUE

!  recover DATAZ from checkpoint file, read year by year
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
! read one line
      READ(30,90,END=1002)IDATE,IYR,NPP(NY,NX) &
        ,(DATAZ(NZ,NY,NX),IsPlantActive(NZ,NY,NX),NZ=1,NPP(NY,NX))
90        FORMAT(2I4,1I3,5(A16,I4))
    ENDDO D9990
  ENDDO D9995

  IF(IDATE.LT.IDAYR.OR.IYR.LT.IYRR)THEN
! mot reaching the recovery date, keep reading
    GO TO 8000
  ELSEIF(IDATE.GE.IDAYR.AND.IYR.EQ.IYRR)THEN
!
!       MATCH PREVIOUS AND CURRENT PLANT SPECIES
! make of copy using pft information from the topographic unit just read in
    D7975: DO NX=NHW,NHE
      D7970: DO NY=NVN,NVS
        D7965: DO NZ=1,NS
          DATAA(NZ,NY,NX)=DATAX(NZ)
          DATAB(NZ,NY,NX)=DATAY(NZ)
        ENDDO D7965
      ENDDO D7970
    ENDDO D7975

    D7995: DO NX=NH1,NH2
      D7990: DO NY=NV1,NV2
        NP(NY,NX)=MAX(NS,NPP(NY,NX))
        D195: DO NN=1,NP(NY,NX)
          DATAP(NN,NY,NX)='NO'
          DATAM(NN,NY,NX)='NO'
        ENDDO D195
! assign the values fo active pft from the record
        IF(NPP(NY,NX).GT.0)THEN
! the check point file has non-zero pft
          D200: DO NN=1,NPP(NY,NX)
            D205: DO NZ=1,NS
              IF(DATAZ(NN,NY,NX).EQ.DATAX(NZ).AND.IsPlantActive(NN,NY,NX).EQ.iPlantIsActive)THEN
                DATAP(NN,NY,NX)=DATAX(NZ)
                DATAM(NN,NY,NX)=DATAY(NZ)
                DATAA(NZ,NY,NX)='NO'
                DATAB(NZ,NY,NX)='NO'
              ENDIF
            ENDDO D205
          ENDDO D200
!
!             ADD NEW PLANT SPECIES
!
          D250: DO NN=1,NP(NY,NX)
            IF(DATAP(NN,NY,NX).EQ.'NO')THEN
              D255: DO NZ=1,NS
                IF(DATAA(NZ,NY,NX).NE.'NO')THEN
                  DATAP(NN,NY,NX)=DATAA(NZ,NY,NX)
                  DATAM(NN,NY,NX)=DATAB(NZ,NY,NX)
                  DATAA(NZ,NY,NX)='NO'
                  DATAB(NZ,NY,NX)='NO'
                  EXIT
                ENDIF
              ENDDO D255
            ENDIF
          ENDDO D250

          D201: DO NZ=NP(NY,NX)+1,JP
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          ENDDO D201
        ELSE

! assign using information from pft files
          D265: DO NZ=1,NS
            DATAP(NZ,NY,NX)=DATAX(NZ)
            DATAM(NZ,NY,NX)=DATAY(NZ)
          ENDDO D265
          D270: DO NZ=NS+1,JP
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          ENDDO D270
        ENDIF
!
!           SET NUMBER OF PLANT SPECIES
!
        NN=JP
        D202: DO NZ=JP,1,-1
          IF(DATAP(NZ,NY,NX).EQ.'NO')THEN
            NN=NN-1
          ELSE
            EXIT
          ENDIF
        ENDDO D202
        NP(NY,NX)=NN
        NP0(NY,NX)=NN
      ENDDO D7990
    ENDDO D7995
  ENDIF
1002 continue
  D6995: DO NX=NHW,NHE
    D6990: DO NY=NVN,NVS
      NP(NY,NX)=NS
      D300: DO NZ=1,NP(NY,NX)
        DATAP(NZ,NY,NX)=DATAX(NZ)
        DATAM(NZ,NY,NX)=DATAY(NZ)
      ENDDO D300
      D301: DO NZ=NP(NY,NX)+1,JP
        DATAP(NZ,NY,NX)='NO'
        DATAM(NZ,NY,NX)='NO'
      ENDDO D301
    ENDDO D6990
  ENDDO D6995

  end subroutine read_checkpt

end module PlantInfoMod
