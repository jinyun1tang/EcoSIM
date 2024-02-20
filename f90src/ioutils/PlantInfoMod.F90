module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use fileUtil, only : open_safe, check_read,int2str,getavu, relavu, opnfil
  use minimathmod, only : isLeap
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use CanopyDataType
  use PlantTraitDataType
  use PlantMathFuncMod
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
  integer :: IDX,IMO,IYR,IDY,ICUT,IDYE,IDYG,IDYS,NumOfCanopyLayersUT,LPY

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
  integer :: NX,NY,NZ,nu_plt
  character(len=128) :: fnm_loc

! begin_execution
  if(disp_planttrait)then
    nu_plt=getavu()
    write(fnm_loc,'(A,I4,A)')'plant_trait.',etimer%get_curr_yearAD(),'.desc'
    call opnfil(fnm_loc,nu_plt,'f')    
  endif
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      D9985: DO NZ=1,NP(NY,NX)
!
!       OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
!       FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!       each column has its own management
!       PREFIX=path for files in current or higher level directory
!       DATAP=PFT file name
        call ReadPlantProperties(nu_plt,NZ,NY,NX)

      ENDDO D9985
    ENDDO D9990
  ENDDO D9995

  if(disp_planttrait)call relavu(nu_plt)
  call ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS)
  RETURN
  END SUBROUTINE readq
!------------------------------------------------------------------------------------------

  subroutine ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS)

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
            iHarvstType_pft(NZ,M,NY,NX)=-1
            jHarvst_pft(NZ,M,NY,NX)=0
            HVST(NZ,M,NY,NX)=1.0E+06_r8
            THIN_pft(NZ,M,NY,NX)=-1.0_r8
            EHVST(1,iplthvst_leaf,NZ,M,NY,NX)=1.0_r8
            EHVST(1,iplthvst_finenonleaf,NZ,M,NY,NX)=1.0_r8
            EHVST(1,iplthvst_woody,NZ,M,NY,NX)=1.0_r8
            EHVST(1,iplthvst_stdead,NZ,M,NY,NX)=1.0_r8
            EHVST(2,iplthvst_leaf,NZ,M,NY,NX)=1.0_r8
            EHVST(2,iplthvst_finenonleaf,NZ,M,NY,NX)=1.0_r8
            EHVST(2,iplthvst_woody,NZ,M,NY,NX)=1.0_r8
            EHVST(2,iplthvst_stdead,NZ,M,NY,NX)=1.0_r8
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
              iDayPlanting_pft(NZ,NY,NX)=IDY
              IYR=yearc
              iYearPlanting_pft(NZ,NY,NX)=MIN(IYR,iYearCurrent)
              IDAYX(NZ,NY,NX)=iDayPlanting_pft(NZ,NY,NX) !planting day
              IYRX(NZ,NY,NX)=iYearPlanting_pft(NZ,NY,NX)   !planting year
              PPZ(NZ,NY,NX)=PPI(NZ,NY,NX)     !population density
            ENDIF

            if(pft_nmgnt(NZ)>0)then
              NN=0
              DO nn1=1,pft_nmgnt(NZ)
                tstr=trim(pft_mgmtinfo(NN1,NZ))
                read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
                READ(TSTR,*)DY,ICUT,NumOfCanopyLayersUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13,&
                   ECUT14,ECUT21,ECUT22,ECUT23,ECUT24

                if(isLeap(iyr) .and. IMO.GT.2)LPY=1
                !obtain the ordinal day
                IF(IMO.EQ.1)then
                  IDY=IDX
                else
                  IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
                endif

                IF(IDY.GT.0.AND.NumOfCanopyLayersUT.EQ.1)THEN
                  iDayPlantHarvest_pft(NZ,NY,NX)=IDY
                  IYR=yearc
                  iYearPlantHarvest_pft(NZ,NY,NX)=MIN(IYR,iYearCurrent)
                ENDIF

                iHarvstType_pft(NZ,IDY,NY,NX)=ICUT
                jHarvst_pft(NZ,IDY,NY,NX)=NumOfCanopyLayersUT
                HVST(NZ,IDY,NY,NX)=HCUT
                THIN_pft(NZ,IDY,NY,NX)=PCUT
                EHVST(1,iplthvst_leaf,NZ,IDY,NY,NX)=ECUT11
                EHVST(1,iplthvst_finenonleaf,NZ,IDY,NY,NX)=ECUT12
                EHVST(1,iplthvst_woody,NZ,IDY,NY,NX)=ECUT13
                EHVST(1,iplthvst_stdead,NZ,IDY,NY,NX)=ECUT14
                EHVST(2,iplthvst_leaf,NZ,IDY,NY,NX)=ECUT21
                EHVST(2,iplthvst_finenonleaf,NZ,IDY,NY,NX)=ECUT22
                EHVST(2,iplthvst_woody,NZ,IDY,NY,NX)=ECUT23
                EHVST(2,iplthvst_stdead,NZ,IDY,NY,NX)=ECUT24

                IF(iHarvstType_pft(NZ,IDY,NY,NX).EQ.4.OR.iHarvstType_pft(NZ,IDY,NY,NX).EQ.6)THEN
                  !animal or insect biomass
                  NN=NN+1
                  if(mod(nn,2)==0)then
                    IDYE=IDY
                    D580: DO IDYG=IDYS+1,IDYE-1
                      iHarvstType_pft(NZ,IDYG,NY,NX)=ICUT
                      jHarvst_pft(NZ,IDYG,NY,NX)=NumOfCanopyLayersUT
                      HVST(NZ,IDYG,NY,NX)=HCUT
                      THIN_pft(NZ,IDYG,NY,NX)=PCUT
                      EHVST(1,iplthvst_leaf,NZ,IDYG,NY,NX)=ECUT11
                      EHVST(1,iplthvst_finenonleaf,NZ,IDYG,NY,NX)=ECUT12
                      EHVST(1,iplthvst_woody,NZ,IDYG,NY,NX)=ECUT13
                      EHVST(1,iplthvst_stdead,NZ,IDYG,NY,NX)=ECUT14
                      EHVST(2,iplthvst_leaf,NZ,IDYG,NY,NX)=ECUT21
                      EHVST(2,iplthvst_finenonleaf,NZ,IDYG,NY,NX)=ECUT22
                      EHVST(2,iplthvst_woody,NZ,IDYG,NY,NX)=ECUT23
                      EHVST(2,iplthvst_stdead,NZ,IDYG,NY,NX)=ECUT24
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
        IDAYY(NZ,NY,NX)=iDayPlantHarvest_pft(NZ,NY,NX)
        IYRY(NZ,NY,NX)=iYearPlantHarvest_pft(NZ,NY,NX)
      ENDDO
    ENDDO
  ENDDO
  end subroutine ReadPlantManagementNC
!------------------------------------------------------------------------------------------

  subroutine ReadPlantProperties(nu_plt,NZ,NY,NX)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: N, NB
  real(r8) :: VRNXI,VRNLI

! begin_execution

!
! READ INPUTS FOR EACH PLANT SPECIES
!
  IF(DATAP(NZ,NY,NX).NE.'NO')THEN
!    write(101,*)'NZ=',NZ,DATAP(NZ,NY,NX)
    call ReadPlantTraitsNC(nu_plt,NZ,NY,NX,VRNLI,VRNXI)
!
!   RE-CALCULATE PLANT INPUTS IN MODEL UNITS
!
!   CanopySWabsorpty_pft,CanopyPARabsorpty_pft=leaf SW,PAR absorbtivity
!
    VmaxRubCarboxyRef_pft(NZ,NY,NX)=2.5_r8*VmaxRubCarboxyRef_pft(NZ,NY,NX)
    VmaxRubOxyRef_pft(NZ,NY,NX)=2.5_r8*VmaxRubOxyRef_pft(NZ,NY,NX)
    VmaxPEPCarboxyRef_pft(NZ,NY,NX)=2.5_r8*VmaxPEPCarboxyRef_pft(NZ,NY,NX)
    SpecChloryfilAct_pft(NZ,NY,NX)=2.5_r8*SpecChloryfilAct_pft(NZ,NY,NX)
    CanopySWabsorpty_pft(NZ,NY,NX)=1.0_r8-CanopySWAlbedo_pft(NZ,NY,NX)-TAUR(NZ,NY,NX)
    CanopyPARabsorpty_pft(NZ,NY,NX)=1.0_r8-CanopyPARalbedo_pft(NZ,NY,NX)-TAUP(NZ,NY,NX)
    CanopySWAlbedo_pft(NZ,NY,NX)=CanopySWAlbedo_pft(NZ,NY,NX)/CanopySWabsorpty_pft(NZ,NY,NX)
    CanopyPARalbedo_pft(NZ,NY,NX)=CanopyPARalbedo_pft(NZ,NY,NX)/CanopyPARabsorpty_pft(NZ,NY,NX)
    TAUR(NZ,NY,NX)=TAUR(NZ,NY,NX)/CanopySWabsorpty_pft(NZ,NY,NX)
    TAUP(NZ,NY,NX)=TAUP(NZ,NY,NX)/CanopyPARabsorpty_pft(NZ,NY,NX)
    SineBranchAngle_pft(NZ,NY,NX)=SIN(BranchAngle_pft(NZ,NY,NX)*RadianPerDegree)
    SinePetioleAngle_pft(NZ,NY,NX)=SIN(PetioleAngle_pft(NZ,NY,NX)*RadianPerDegree)
    MatureGroup_pft(NZ,NY,NX)=GROUPX(NZ,NY,NX)

    IF(iPlantTurnoverPattern_pft(NZ,NY,NX).NE.0)THEN
!
      RefNodeInitRate_pft(NZ,NY,NX)=RefNodeInitRate_pft(NZ,NY,NX)/MaxNodesPerBranch
      RefLeafAppearRate_pft(NZ,NY,NX)=RefLeafAppearRate_pft(NZ,NY,NX)/MaxNodesPerBranch
      MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)/MaxNodesPerBranch
      XTLI(NZ,NY,NX)=XTLI(NZ,NY,NX)/MaxNodesPerBranch
    ENDIF
    MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)-XTLI(NZ,NY,NX)
    IF(CriticalPhotoPeriod_pft(NZ,NY,NX).LT.0.0_r8)THEN
      CriticalPhotoPeriod_pft(NZ,NY,NX)=DayLenthMax(NY,NX)
    ENDIF
    D5: DO NB=1,NumOfCanopyLayers
      IF(iPlantPhenologyType_pft(NZ,NY,NX).EQ.iphenotyp_evgreen.AND.iPlantPhenologyPattern_pft(NZ,NY,NX).NE.iplt_annual)THEN
        HourThreshold4LeafOut_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNLI+144.0_r8*iPlantInitThermoAdaptZone(NZ,NY,NX)*(NB-1))
        HourThreshold4LeafOff_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNXI+144.0_r8*iPlantInitThermoAdaptZone(NZ,NY,NX)*(NB-1))
      ELSE
        HourThreshold4LeafOut_brch(NB,NZ,NY,NX)=VRNLI
        HourThreshold4LeafOff_brch(NB,NZ,NY,NX)=VRNXI
      ENDIF
    ENDDO D5
  ENDIF
! WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,CriticalPhotoPeriod_pft(NZ,NY,NX)
!1111    FORMAT(A20,2I8,E12.4)
  end subroutine ReadPlantProperties

!------------------------------------------------------------------------------------------
  subroutine ReadPlantTraitsNC(nu_plt,NZ,NY,NX,VRNLI,VRNXI)

  use GrosubPars, only : get_pft_loc
  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  real(r8), intent(out) :: VRNLI,VRNXI
  character(len=40) :: pft_lname
  character(len=64):: koppen_climl
  character(len=3) :: koppen_clims
  integer :: loc,N

  loc=get_pft_loc(DATAP(NZ,NY,NX)(1:6),pft_lname,koppen_climl,koppen_clims)
  DATAPI(NZ,NY,NX)=loc
  call ncd_getvar(pft_nfid, 'ICTYP', loc, iPlantPhotosynthesisType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IGTYP', loc, iPlantMorphologyType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ISTYP', loc, iPlantPhenologyPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IDTYP', loc, iPlantDevelopPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'INTYP', loc, iPlantNfixType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IWTYP', loc, iPlantPhenologyType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IPTYP', loc, iPlantPhotoperiodType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IBTYP', loc, iPlantTurnoverPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IRTYP', loc, iPlantGrainType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'MY', loc, MY(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ZTYPI', loc, iPlantInitThermoAdaptZone(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'VCMX', loc, VmaxRubCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VOMX', loc, VmaxRubOxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VCMX4', loc, VmaxPEPCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO2', loc, XKCO2(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKO2', loc, XKO2(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO24', loc, Km4PEPCarboxy_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RUBP', loc, LeafRuBPConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PEPC', loc, FracLeafProtinAsPEPCarboxyl_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ETMX', loc, SpecChloryfilAct_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL', loc, LeafC3ChlorofilConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL4', loc, LeafC4ChlorofilConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'FCO2', loc, CanPCi2CaRatio(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'ALBR', loc, CanopySWAlbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ALBP', loc, CanopyPARalbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUR', loc, TAUR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUP', loc, TAUP(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'XRNI', loc, RefNodeInitRate_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XRLA', loc, RefLeafAppearRate_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CTC', loc, TCelciusChill4Seed(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VRNLI', loc,VRNLI)
  call ncd_getvar(pft_nfid, 'VRNXI', loc,VRNXI)
  call ncd_getvar(pft_nfid, 'WDLF', loc,WDLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PB', loc,MinNonstructalC4InitBranch(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'GROUPX', loc,GROUPX(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XTLI', loc,XTLI(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XDL', loc,CriticalPhotoPeriod_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XPPD', loc,PhotoPeriodSens_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'SLA1', loc,SLA1(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SSL1', loc,PetoLen2Mass_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SNL1', loc,SNL1(NZ,NY,NX))


  call ncd_getvar(pft_nfid, 'CLASS', loc,CLASS(1:NumOfLeafZenithSectors,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CFI', loc,ClumpFactorInit_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGBR', loc,BranchAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGSH', loc,PetioleAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'STMX', loc,MaxPotentSeedNumber_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SDMX', loc,MaxSeedNumPerSite_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRMX', loc,MaxSeedCMass(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRDM', loc,SeedCMass(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GFILL', loc,GrainFillRateat25C_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'WTSTDI', loc,StandingDeadInitC_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'RRAD1M', loc,Max1stRootRadius(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RRAD2M', loc,Max2ndRootRadius(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PORT', loc,RootPorosity(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PR', loc,MinNonstructuralC4InitRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRR', loc,RSRR(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRA', loc,RSRA(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PTSHT', loc,ShutRutNonstructElmntConducts_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RTFQ', loc,RootBranchFreq_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZH', loc,VmaxNH4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZH', loc,KmNH4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZH', loc,CMinNH4Root_pft(ipltroot,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZO', loc,VmaxNO3Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZO', loc,KmNO3Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZO', loc,CminNO3Root_pft(ipltroot,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXPO', loc,VmaxPO4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMPO', loc,KmPO4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNPO', loc,CMinPO4Root_pft(ipltroot,NZ,NY,NX))

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
  call ncd_getvar(pft_nfid, 'DMND', loc,NoduGrowthYield_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CNLF', loc,CNLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSHE', loc,CNSHE(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSTK', loc,rNCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRSV', loc,rNCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNHSK', loc,rNCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNEAR', loc,rNCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNGR', loc,CNGR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRT', loc,RootrNC_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNND', loc,NodulerNC_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CPLF', loc,CPLF(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSHE', loc,CPSHE(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSTK', loc,rPCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRSV', loc,rPCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPHSK', loc,rPCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPEAR', loc,rPCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPGR', loc,CPGR(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRT', loc,RootrPC_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPND', loc,NodulerPC_pft(NZ,NY,NX))

  if(disp_planttrait)then
    call pft_display(nu_plt,NZ,NY,NX,pft_lname,koppen_climl,koppen_clims)
    call photosyns_trait_disp(nu_plt,NZ,NY,NX)
    call plant_optic_trait_disp(nu_plt,NZ,NY,NX)
    call Phenology_trait_disp(nu_plt,NZ,NY,NX,VRNLI,VRNXI)
    call morphology_trait_disp(nu_plt,NZ,NY,NX)
    call Root_trait_disp(nu_plt,NZ,NY,NX)
    call Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
    call plant_water_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  endif

  end subroutine ReadPlantTraitsNC
!------------------------------------------------------------------------------------------
  subroutine pft_display(nu_plt,NZ,NY,NX,pft_lname,koppen_climl,koppen_clims)
  use abortutils , only : endrun  
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  character(len=*),intent(in) :: pft_lname  
  character(len=*),intent(in) :: koppen_climl
  character(len=*),intent(in) :: koppen_clims
  integer :: j

!   iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4
!   iPlantMorphologyType_pft=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   iPlantPhenologyPattern_pft=growth habit:0=annual,1=perennial
!   iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate
!   iPlantNfixType=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   iPlantPhenologyType_pft=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   iPlantPhotoperiodType_pft=photoperiod type:0=day neutral,1=short day,2=long day
!   iPlantTurnoverPattern_pft=turnover:if iPlantMorphologyType_pft=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if iPlantMorphologyType_pft=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   iPlantGrainType_pft=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   iPlantInitThermoAdaptZone=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical
  write(nu_plt,*)('=',j=1,100)
  write(nu_plt,*)'PLANT traits for FUNCTIONAL TYPE (NZ,NY,NX)=',NZ,NY,NX,DATAP(NZ,NY,NX)(1:6)
  write(nu_plt,*)'Plant name: ',pft_lname
  write(nu_plt,*)'koppen climate info:',koppen_clims//','//koppen_climl
  select CASE (iPlantPhotosynthesisType(NZ,NY,NX))
  case (3)
    write(nu_plt,*)'C3 photosynthesis'
  case (4)
    write(nu_plt,*)'C4 photosynthesis'
  case default
    write(nu_plt,*)'photosynthesis type not defined'
  end select

  select case(iPlantMorphologyType_pft(NZ,NY,NX))
  case (0)
    write(nu_plt,*)'shallow root profile, like bryophytes'
  case (1)
    write(nu_plt,*)'intermediate root profile, like herbs'
  case (2)
    write(nu_plt,*)'deep root profile, like trees'
  case default
    write(nu_plt,*)'root profile not defined'
  end select

  select case (iPlantPhenologyPattern_pft(NZ,NY,NX))
  case (0)
    write(nu_plt,*)'Annual plant'
  case (1)
    write(nu_plt,*)'perennial plant'
  case default
    write(nu_plt,*)'growth habit not defined'
  end select

  select case (iPlantDevelopPattern_pft(NZ,NY,NX))
  case (0)
    write(nu_plt,*)'determinate growth pattern'
  case (1)
    write(nu_plt,*)'indetermimate growth pattern'
  case default
    write(nu_plt,*)'growth pattern not defined'
  end select

  select case (iPlantNfixType(NZ,NY,NX))
! 1,2, 3, e.g. legumes
  case (in2fixtyp_root_fast)
    write(nu_plt,*)'Rapid root N2 fixation symbiosis'
  case (in2fixtyp_root_medium)
    write(nu_plt,*)'Intermediate root N2 fixation symbiosis'
  case (in2fixtyp_root_slow)
    write(nu_plt,*)'Slow root N2 fixation symbiosis'
!4,5,6, e.g. cyanobacteria
  case (in2fixtyp_canopy_fast)
    write(nu_plt,*)'Rapid canopy N2 fixation symbiosis'
  case (in2fixtyp_canopy_medium)
    write(nu_plt,*)'Intermediate canopy N2 fixation symbiosis'
  case (in2fixtyp_canopy_slow)
    write(nu_plt,*)'Slow canopy N2 fixation symbiosis'
  case default
    write(nu_plt,*)'No N2 fixation symbiosis defined'
  end select

  select case(iPlantPhenologyType_pft(NZ,NY,NX))
  case (iphenotyp_evgreen)
    write(nu_plt,*)'phenology type: evergreen'
  case (iphenotyp_coldecidu)
    write(nu_plt,*)'phenology type: cold deciduous'
  case (iphenotyp_drouhtdecidu)
    write(nu_plt,*)'phenology type: drought deciduous'
  case (iphenotyp_coldroutdecidu)
    write(nu_plt,*)'phenology type: cold+drought deciduous'
  case default
    write(nu_plt,*)'phenology type not defined'
  end select

  select case(iPlantPhotoperiodType_pft(NZ,NY,NX))
  case (iphotop_neutral)
    write(nu_plt,*)'day neutral photoperiod'
  case (iphotop_short)
    write(nu_plt,*)'short day photoperiod'
  case (iphotop_long)
    write(nu_plt,*)'long day photoperiod'
  case default
    write(nu_plt,*)'photoperiod not defined'
  end select

  if(is_plant_treelike(iPlantMorphologyType_pft(NZ,NY,NX)))then
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(nu_plt,*)'Rapid plant biome turnover (deciduous)'
    case (2)
      write(nu_plt,*)'Very slow plant biome turnover (needleleaf evergreen)'
    case (3)
      write(nu_plt,*)'Slow plant biome turnover (Broadleaf evergreen)'
    case (4)
      write(nu_plt,*)'Plant biome turnover semi-deciduous'
    case (5)
      write(nu_plt,*)'Plant biome turnover semi-evergreen'
    case default
      write(nu_plt,*)'Plant biome turnover not defined'
    end select
  else
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(nu_plt,*)'Rapid all aboveground plant biome turnover (herbaceous)'
    case default
      write(nu_plt,*)'Plant biome turnover not defined'
    end select
  endif

  select case(iPlantGrainType_pft(NZ,NY,NX))
  case (igraintyp_abvgrnd)
    write(nu_plt,*)'Above ground storage organ'
  case (igraintyp_blwgrnd)
    write(nu_plt,*)'Belowground storage organ'
  case default
    write(nu_plt,*)'Storage organ not defined'
  end select

  select case(MY(NZ,NY,NX))
  case (1)
    write(nu_plt,*)'No mycorrhizal'
  case (2)
    write(nu_plt,*)'Mycorrhizal'
  case default
    write(nu_plt,*)'Wrong option for mycorrhizae'
  end select

  select case(INT(iPlantInitThermoAdaptZone(NZ,NY,NX)))
  case (ithermozone_arcboreal)
    write(nu_plt,*)'thermal adaptation zone: arctic, boreal'
  case (ithermozone_cooltempr)
    write(nu_plt,*)'thermal adaptation zone: cool temperate'
  case (ithermozone_warmtempr)
    write(nu_plt,*)'thermal adaptation zone: warm temperate'
  case (ithermozone_subtropic)
    write(nu_plt,*)'thermal adaptation zone: subtropical'
  case (ithermozone_tropical)
    write(nu_plt,*)'thermal adaptation zone: tropical'
  case default
    write(nu_plt,*)'Not thermal adaptation zone defined'
  end select

  end subroutine pft_display

!------------------------------------------------------------------------------------------
  subroutine photosyns_trait_disp(nu_plt,NZ,NY,NX)
  !
  !Description
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'PHOTOSYNTHETIC PROPERTIES'
  write(nu_plt,*)'Specific rubisco carboxylase (umol C g-1 s-1): ',VmaxRubCarboxyRef_pft(NZ,NY,NX)
  write(nu_plt,*)'Specific rubisco oxygenase (umol O2 g-1 s-1): ',VmaxRubOxyRef_pft(NZ,NY,NX)
  write(nu_plt,*)'Specific PEP carboxylase activity (umol g-1 s-1): '&
    ,VmaxPEPCarboxyRef_pft(NZ,NY,NX)
  write(nu_plt,*)'Km for VmaxRubCarboxyRef_pft(uM): ',XKCO2(NZ,NY,NX)
  write(nu_plt,*)'Km for VmaxRubOxyRef_pft (uM): ',XKO2(NZ,NY,NX)
  write(nu_plt,*)'KM for VmaxPEPCarboxyRef_pft (uM): ',Km4PEPCarboxy_pft(NZ,NY,NX)
  write(nu_plt,*)'Fraction of leaf protein in rubisco (g rub/(g protein)): ',LeafRuBPConc_pft(NZ,NY,NX)
  write(nu_plt,*)'Fraction of leaf protein in PEP carboxylase (g pep/(g protein)): ',FracLeafProtinAsPEPCarboxyl_pft(NZ,NY,NX)
  write(nu_plt,*)'Specific chlorophyll activity (umol e- gC-1 s-1): ',SpecChloryfilAct_pft(NZ,NY,NX)
  if(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic3_photo)then
    write(nu_plt,*)'Fraction of leaf protein as chlorophyll in mesophyll (C3) (g Chl /(g protein))',LeafC3ChlorofilConc_pft(NZ,NY,NX)
  elseif(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic4_photo)then
    write(nu_plt,*)'Fraction of leaf protein as chlorophyll in bundle sheath(C4) (g Chl /(g protein)):',LeafC3ChlorofilConc_pft(NZ,NY,NX)
    write(nu_plt,*)'fraction of leaf protein in mesophyll chlorophyll(C4) (g Chl /(g protein)):',LeafC4ChlorofilConc_pft(NZ,NY,NX)
  endif
  write(nu_plt,*)'intercellular:atmospheric CO2 concentration ratio:',CanPCi2CaRatio(NZ,NY,NX)
  end subroutine photosyns_trait_disp

!------------------------------------------------------------------------------------------

  subroutine Phenology_trait_disp(nu_plt,NZ,NY,NX,VRNLI,VRNXI)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: VRNLI,VRNXI
  integer, intent(in) :: nu_plt
  integer :: J

  write(nu_plt,*)('-',j=1,100)

  write(nu_plt,*)'PHENOLOGICAL PROPERTIES'
  write(nu_plt,*)'rate of node initiation at 25oC (h-1): XRNI',RefNodeInitRate_pft(NZ,NY,NX)
  write(nu_plt,*)'rate of leaf appearance at 25oC (h-1): XRLA',RefLeafAppearRate_pft(NZ,NY,NX)
  write(nu_plt,*)'chilling temperature for CO2 fixation, '// &
    'seed loss (oC): CTC',TCelciusChill4Seed(NZ,NY,NX)
  write(nu_plt,*)'hour requirement for spring leafout: VRNLI',VRNLI
  write(nu_plt,*)'hour requirement for autumn leafoff: VRNXI',VRNXI
  write(nu_plt,*)'leaf length:width ratio: WDLF',WDLF(NZ,NY,NX)
  write(nu_plt,*)'nonstructural C concentration needed for branching'// &
    ':PB',MinNonstructalC4InitBranch(NZ,NY,NX)
  end subroutine Phenology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine morphology_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: N,j

  write(nu_plt,*)('-',j=1,100)  
  write(nu_plt,*)'MORPHOLOGICAL PROPERTIES'
  write(nu_plt,*)'growth in leaf area vs mass (m2 g-1): SLA1 ',SLA1(NZ,NY,NX)
  write(nu_plt,*)'growth in petiole length vs mass (m g-1): SSL1 ',PetoLen2Mass_pft(NZ,NY,NX)
  write(nu_plt,*)'growth in internode length vs mass (m g-1): SNL1',SNL1(NZ,NY,NX)
  write(nu_plt,*)'fraction of leaf area in 0-22.5,45,67.5,90o '// &
    'inclination classes: CLASS',(CLASS(N,NZ,NY,NX),N=1,NumOfLeafZenithSectors)
  write(nu_plt,*)'initial clumping factor: CFI',ClumpFactorInit_pft(NZ,NY,NX)
  write(nu_plt,*)'stem angle from horizontal: BranchAngle_pft',BranchAngle_pft(NZ,NY,NX)
  write(nu_plt,*)'petiole angle from horizontal: PetioleAngle_pft',PetioleAngle_pft(NZ,NY,NX)
  write(nu_plt,*)'maximum potential seed mumber from '// &
    'pre-anthesis stalk growth: STMX',MaxPotentSeedNumber_pft(NZ,NY,NX)
  write(nu_plt,*)'maximum seed number per MaxPotentSeedNumber_pft: SDMX',MaxSeedNumPerSite_pft(NZ,NY,NX)
  write(nu_plt,*)'maximum seed size per MaxSeedNumPerSite_pft (g): GRMX',MaxSeedCMass(NZ,NY,NX)
  write(nu_plt,*)'seed size at planting (gC): GRDM',SeedCMass(NZ,NY,NX)    !could be greater than MaxSeedCMass, accouting for seedling
  write(nu_plt,*)'grain filling rate at 25 oC (g seed-1 h-1): GFILL',GrainFillRateat25C_pft(NZ,NY,NX)
  write(nu_plt,*)'mass of dead standing biomass at planting (gC m-2): WTSTDI',StandingDeadInitC_pft(NZ,NY,NX)
  end subroutine morphology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ROOT CHARACTERISTICS'
  write(nu_plt,*)'radius of primary roots (m): Max1stRootRadius',Max1stRootRadius(1,NZ,NY,NX)
  write(nu_plt,*)'radius of secondary roots (m): Max2ndRootRadius',Max2ndRootRadius(1,NZ,NY,NX)
  write(nu_plt,*)'primary/fine root porosity (m3 m-3): PORT',RootPorosity(1,NZ,NY,NX)
  write(nu_plt,*)'nonstructural C concentration needed for root'// &
    ' branching (gC/gC): PR',MinNonstructuralC4InitRoot_pft(NZ,NY,NX)
  write(nu_plt,*)'radial root resistivity for water uptake (m2 MPa-1 h-1): RSRR',RSRR(1,NZ,NY,NX)
  write(nu_plt,*)'axial root resistivity for water uptake (m2 MPa-1 h-1): RSRA',RSRA(1,NZ,NY,NX)
  write(nu_plt,*)'rate constant for equilibrating shoot-root '// &
    'nonstructural C concn: PTSHT',ShutRutNonstructElmntConducts_pft(NZ,NY,NX)
  write(nu_plt,*)'root branching frequency (m-1): RootBranchFreq_pft',RootBranchFreq_pft(NZ,NY,NX)
  end subroutine Root_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ROOT UPTAKE PARAMETERS'
  write(nu_plt,*)'NH4 max uptake (g m-2 h-1): UPMXZH',VmaxNH4Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'NH4 uptake Km (uM): UPKMZH',KmNH4Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'NH4 uptake min concn (uM): UPMNZH',CMinNH4Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'NO3 max uptake (g m-2 h-1): UPMXZO',VmaxNO3Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'NO3 uptake Km (uM): UPKMZO',KmNO3Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'NO3 uptake min concn (uM): UPMNZO',CminNO3Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'H2PO4 or H1PO4 max uptake (g m-2 h-1): UPMXPO',VmaxPO4Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'H2PO4 or H1PO4 uptake Km (uM): UPKMPO',KmPO4Root_pft(ipltroot,NZ,NY,NX)
  write(nu_plt,*)'H2PO4 or H1PO4 uptake min conc (uM): UPMNPO',CMinPO4Root_pft(ipltroot,NZ,NY,NX)
  end subroutine Root_nutrient_trait_disp


!------------------------------------------------------------------------------------------
  subroutine plant_water_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)  
  write(nu_plt,*)'WATER RELATIONS'
  write(nu_plt,*)'leaf osmotic potential at zero leaf water '// &
     'potential (MPa): OSMO',OSMO(NZ,NY,NX)
  write(nu_plt,*)'shape parameter for stomatal resistance vs '// &
     'leaf turgor potential: RCS',RCS(NZ,NY,NX)
  write(nu_plt,*)'cuticular resistance (s m-1): RSMX',RSMX(NZ,NY,NX)
  end subroutine plant_water_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ORGAN GROWTH YIELDS'
  write(nu_plt,*)'leaf dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMLF',LeafBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'petiole dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMSHE',PetioleBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'stalk dry matter C production vs '// &
    'nonstructural C consumption (g g-1): DMSTK',StalkBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'stalk reserve C production vs '// &
    'nonstructural C consumption (g g-1): DMRSV',ReserveBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'husk dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1): DMHSK',HuskBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'ear dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1): DMEAR',EarBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'grain C production vs nonstructural C'// &
    ' consumption (g g-1): DMGR',GrainBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'root dry matter C production vs nonstructural C'// &
    ' consumption (g g-1): DMRT',RootBiomGrowthYield(NZ,NY,NX)
  write(nu_plt,*)'nodule bacteria in root nodule,canopy dry matter'// &
    'C production vs nonstructural C consumption (g g-1): DMND' &
    ,NoduGrowthYield_pft(NZ,NY,NX)
  end subroutine plant_biomyield_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ORGAN N AND P CONCENTRATIONS'
  write(nu_plt,*)'NC ratio in plant leaves (gN/gC): CNLF',CNLF(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant petiole (gN/gC): CNSHE',CNSHE(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant stalk (gN/gC): CNSTK',rNCStalk_pft(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant stalk reserve (gN/gC): CNRSV',rNCReserve_pft(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant husk (gN/gC): CNHSK',rNCHusk_pft(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant ear (gN/gC): CNEAR',rNCEar_pft(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant grain (gN/gC): CNGR',CNGR(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant root (gN/gC): CNumRootAxes_pft',RootrNC_pft(NZ,NY,NX)
  write(nu_plt,*)'NC ratio in plant nodule (gN/gC): CNND',NodulerNC_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant leaves (gP/gC): CPLF',CPLF(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant petiole (gP/gC): CPSHE',CPSHE(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant stalk (gP/gC): CPSTK',rPCStalk_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant stalk reserve (gP/gC): CPRSV',rPCReserve_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant husk (gP/gC): CPHSK',rPCHusk_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant ear (gP/gC): CPEAR',rPCEar_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant grain (gP/gC): CPGR',CPGR(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant root (gP/gC): CPRT',RootrPC_pft(NZ,NY,NX)
  write(nu_plt,*)'PC ratio in plant nodule (gP/gC): CPND',NodulerPC_pft(NZ,NY,NX)
  end subroutine plant_biomstoich_trait_disp

!------------------------------------------------------------------------------------------

  subroutine plant_optic_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'OPTICAL PROPERTIES'
  write(nu_plt,*)'leaf SW albedo: ALBR',CanopySWAlbedo_pft(NZ,NY,NX)
  write(nu_plt,*)'leaf PAR albedo: ALBP',CanopyPARalbedo_pft(NZ,NY,NX)
  write(nu_plt,*)'leaf SW transmission: TAUR',TAUR(NZ,NY,NX)
  write(nu_plt,*)'leaf PAR transmission: TAUP',TAUP(NZ,NY,NX)
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
          iDayPlanting_pft(NZ,NY,NX)=-1E+06
          iYearPlanting_pft(NZ,NY,NX)=-1E+06
          iDayPlantHarvest_pft(NZ,NY,NX)=1E+06
          iYearPlantHarvest_pft(NZ,NY,NX)=1E+06
        enddo
      enddo
    ENDDO D9999

  ENDIF
  if(lverb)write(*,*)'ReadPlantInfoNC'
  call ReadPlantInfoNC(yeari,NE,NEX,NHW,NHE,NVN,NVS)

  END subroutine routq


!------------------------------------------------------------------------------------------

  subroutine ReadPlantInfoNC(yeari,NE,NEX,NHW,NHE,NVN,NVS)
  !
  !Description
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
                IF(KoppenClimZone(NY,NX).GT.0)THEN
                  WRITE(CLIMATE,'(I2)')KoppenClimZone(NY,NX)
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
        ,(DATAZ(NZ,NY,NX),IsPlantActive_pft(NZ,NY,NX),NZ=1,NPP(NY,NX))
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
              IF(DATAZ(NN,NY,NX).EQ.DATAX(NZ).AND.IsPlantActive_pft(NN,NY,NX).EQ.iPlantIsActive)THEN
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
