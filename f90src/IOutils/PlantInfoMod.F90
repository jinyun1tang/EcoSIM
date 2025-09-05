module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use fileUtil, only : open_safe, check_read,int2str,getavu, relavu, opnfil,iulog
  use minimathmod, only : isLeap
  use abortutils, only : endrun  
  use PlantTraitTableMod  
  use netcdf
  use ncdio_pio  
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use CanopyDataType
  use PlantTraitDataType
  use PlantMathFuncMod
  use PlantMgmtDataType
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

  public :: ReadPlantInfo
  public :: ReadPlantTraitTable
  contains

  subroutine ReadPlantInfo(yearc,yeari,NHW,NHE,NVN,NVS)
!
! DESCRIPTION
!

  implicit none
  integer, intent(in) :: yearc, yeari
  integer, intent(in) :: NHW,NHE,NVN,NVS     !simulation domain specification

! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
  if(lverb)WRITE(*,333)'ROUTQ'
  CALL ReadPlantInfoNC(yeari,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  if(lverb)WRITE(*,333)'READQ'
  CALL READQ(yearc,yeari,NHW,NHE,NVN,NVS)

333   FORMAT(A8)

  end subroutine ReadPlantInfo
!------------------------------------------------------------------------------------------

  SUBROUTINE readq(yearc,yeari,NHW,NHE,NVN,NVS)
!!
! Description
! THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
!     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
!
  implicit none
  integer, intent(in) :: yearc   !current model year
  integer, intent(in) :: yeari   !current forcing year 
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NX,NY,NZ,nu_plt
  character(len=128) :: fnm_loc
  logical :: pft_changed

! begin_execution
  call ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS,pft_changed)

  if(disp_planttrait .and. pft_changed)then
    nu_plt=getavu()
    write(fnm_loc,'(A,I4,A)')'plant_trait.',etimer%get_curr_yearAD(),'.desc'
    call opnfil(fnm_loc,nu_plt,'f')    
  endif
  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      D9985: DO NZ=1,NP_col(NY,NX)
!
!       OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
!       FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!       each column has its own management
!       PREFIX=path for files in current or higher level directory
!       DATAP=PFT file name
        call ReadPlantProperties(nu_plt,NZ,NY,NX,pft_changed)

      ENDDO D9985
    ENDDO D9990
  ENDDO D9995

  if(disp_planttrait .and. pft_changed)call relavu(nu_plt)

  RETURN
  END SUBROUTINE readq
!------------------------------------------------------------------------------------------
  subroutine readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  !
  !read planting information
  implicit none
  type(file_desc_t), intent(in) :: pftinfo_nfid
  integer, intent(in) :: ntopou
  integer, intent(in) :: iyear    !file record number
  integer, intent(in) :: yearc    !current model year
  integer, intent(in) :: NHW,NHE,NVN,NVS
  logical  :: readvar
  integer  :: NH1,NV1,NH2,NV2,NS
  integer :: NTOPO,NY,NX,NZ
  real(r8) :: DY
  type(Var_desc_t) :: vardesc
  character(len=128) :: pft_pltinfo(JP),tstr
  integer :: LPY,IDX,IMO,IYR,IDY

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        PlantinDepz_pft(NZ,NY,NX)=ppmc
      ENDDO
    ENDDO
  ENDDO

  DO NTOPO=1,ntopou
    call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
    call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)
    if(any((/NH1, NV1,NH2,NV2/)<=0) .or. any((/NH1, NV1,NH2,NV2/)>ntopou))exit
    !planting information
    call check_var(pftinfo_nfid, 'pft_pltinfo', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find pft_pltinfo in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_pltinfo, &
      start = (/1,1,NTOPO,iyear/),count = (/len(pft_pltinfo(1)),JP,1,1/)),&
      trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

    DO NX=NH1,NH2
      DO NY=NV1,NV2
        DO NZ=1,MIN(NS,NP_col(NY,NX))
          tstr=trim(pft_pltinfo(NZ))
          if (tstr .EQ. "") then
            cycle
          endif
          read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
          read(tstr,*)DY,PPI_pft(NZ,NY,NX),PlantinDepz_pft(NZ,NY,NX)
          PlantinDepz_pft(NZ,NY,NX)=AZMAX1(PlantinDepz_pft(NZ,NY,NX),ppmc)
          LPY=0
          if(isLeap(iyr) .and. IMO.GT.2)LPY=1
          !obtain the ordinal day
          IF(IMO.EQ.1)then
            IDY=IDX
          else
            IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
          endif
          IF(IDY.GT.0.AND.IYR.GT.0)THEN
            iDayPlanting_pft(NZ,NY,NX)  = IDY
            IYR                         = yearc
            iYearPlanting_pft(NZ,NY,NX) = MIN(IYR,iYearCurrent)
            iPlantingDay_pft(NZ,NY,NX)  = iDayPlanting_pft(NZ,NY,NX) !planting day
            iPlantingYear_pft(NZ,NY,NX) = iYearPlanting_pft(NZ,NY,NX)   !planting year
            PPatSeeding_pft(NZ,NY,NX)   = PPI_pft(NZ,NY,NX)     !population density
          ENDIF          
        ENDDO
      ENDDO
    ENDDO
  ENDDO  
  end subroutine readplantinginfo  

!------------------------------------------------------------------------------------------
  subroutine InitPlantMgmnt(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,NZ,M

! initialize the disturbance arrays
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        DO M=1,366
          iHarvstType_pft(NZ,M,NY,NX)                                   = -1
          jHarvstType_pft(NZ,M,NY,NX)                                       = 0
          CanopyHeightCut_pft(NZ,M,NY,NX)                               = 1.0E+06_r8
          THIN_pft(NZ,M,NY,NX)                                          = -1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,M,NY,NX)        = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,M,NY,NX) = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_woody,NZ,M,NY,NX)       = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,M,NY,NX)      = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,M,NY,NX)        = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,M,NY,NX) = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ,M,NY,NX)       = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,M,NY,NX)      = 1.0_r8
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine InitPlantMgmnt
!------------------------------------------------------------------------------------------
  subroutine readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  !
  !read plant management information
  implicit none
  type(file_desc_t), intent(in) :: pftinfo_nfid
  integer, intent(in) :: ntopou
  integer, intent(in) :: iyear    !file record number
  integer, intent(in) :: yearc    !current model year  
  integer, intent(in) :: NHW,NHE,NVN,NVS

  logical  :: readvar
  real(r8) :: DY,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
  real(r8) :: ECUT24,HCUT,PCUT
  integer  :: pft_nmgnt(JP)
  integer  :: NH1,NV1,NH2,NV2,NS
  type(Var_desc_t) :: vardesc
  integer :: NTOPO,NY,NX,NZ
  character(len=128) :: pft_pltinfo(JP),tstr  
  character(len=128) :: pft_mgmtinfo(24,JP)
  integer :: LPY,IDX,IMO,IYR,IDY,ICUT,IDYE,IDYG,IDYS
  integer :: M,NN,N,nn1
  integer :: JCUT

  call InitPlantMgmnt(NHW,NHE,NVN,NVS)
  
  DO NTOPO=1,ntopou
    call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
    call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)

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
    if(any((/NH1, NV1,NH2,NV2/)<=0) .or. any((/NH1, NV1,NH2,NV2/)>ntopou))exit

    DO NX=NH1,NH2
      DO NY=NV1,NV2
        DO NZ=1,MIN(NS,NP_col(NY,NX))
          if(pft_nmgnt(NZ)>0)then
            NN=0
            DO nn1=1,pft_nmgnt(NZ)
              if(len_trim(pft_mgmtinfo(NN1,NZ))==0)cycle
              tstr=trim(pft_mgmtinfo(NN1,NZ))
              !read day/month/year when management occurs
              read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
              READ(TSTR,*)DY,ICUT,JCUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13,&
                  ECUT14,ECUT21,ECUT22,ECUT23,ECUT24
              LPY=0
              if(isLeap(iyr) .and. IMO.GT.2)LPY=1
              !obtain the ordinal day
              IF(IMO.EQ.1)then
                IDY=IDX
              else
                IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
              endif

              IF(IDY.GT.0 .AND. JCUT.EQ.1)THEN
                iDayPlantHarvest_pft(NZ,NY,NX)=IDY
                IYR=yearc
                iYearPlantHarvest_pft(NZ,NY,NX)=MIN(IYR,iYearCurrent)
              ENDIF
              !harvest is specified with two type numbers, first is large categroy, second is pft-specific operation
              iHarvstType_pft(NZ,IDY,NY,NX)                         = ICUT
              jHarvstType_pft(NZ,IDY,NY,NX)                         = JCUT
              CanopyHeightCut_pft(NZ,IDY,NY,NX)                     = HCUT
              THIN_pft(NZ,IDY,NY,NX)                                = PCUT
              FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,IDY,NY,NX)        = ECUT11
              FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,IDY,NY,NX) = ECUT12
              FracBiomHarvsted(iHarvst_pft,iplthvst_woody,NZ,IDY,NY,NX)       = ECUT13
              FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,IDY,NY,NX)      = ECUT14
              
              FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,IDY,NY,NX)        = ECUT21
              FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,IDY,NY,NX) = ECUT22
              FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ,IDY,NY,NX)       = ECUT23
              FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,IDY,NY,NX)      = ECUT24

              IF(iHarvstType_pft(NZ,IDY,NY,NX).EQ.4.OR.iHarvstType_pft(NZ,IDY,NY,NX).EQ.6)THEN
                !animal or insect biomass
                NN=NN+1
                if(mod(nn,2)==0)then
                  IDYE=IDY
                  D580: DO IDYG=IDYS+1,IDYE-1
                    iHarvstType_pft(NZ,IDYG,NY,NX)                         = ICUT
                    jHarvstType_pft(NZ,IDYG,NY,NX)                             = JCUT
                    CanopyHeightCut_pft(NZ,IDYG,NY,NX)                     = HCUT
                    THIN_pft(NZ,IDYG,NY,NX)                                = PCUT
                    FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,IDYG,NY,NX)        = ECUT11
                    FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,IDYG,NY,NX) = ECUT12
                    FracBiomHarvsted(iHarvst_pft,iplthvst_woody,NZ,IDYG,NY,NX)       = ECUT13
                    FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,IDYG,NY,NX)      = ECUT14
                    FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,IDYG,NY,NX)        = ECUT21
                    FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,IDYG,NY,NX) = ECUT22
                    FracBiomHarvsted(iHarvst_col,iplthvst_woody,NZ,IDYG,NY,NX)       = ECUT23
                    FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,IDYG,NY,NX)      = ECUT24
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
  end subroutine readplantmgmtinfo
!------------------------------------------------------------------------------------------

  subroutine ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS,pft_changed)
  !
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: yeari   !forcing data year
  integer, intent(in) :: yearc   !current model year
  logical, optional, intent(out):: pft_changed
  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear     !data record year
  integer :: year      !year of the data
  integer :: ntopou,NY,NX,NZ
  logical :: pft_changed_loc

  pft_changed_loc=.false.
  if (len_trim(pft_mgmt_in)==0)return
  
  call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

  call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)

  call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
    trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))

  !obtain the number of topographic units
  ntopou=get_dim_len(pftinfo_nfid, 'ntopou')    
  if(first_topou)ntopou=1
  
  if(pft_dflag==0)then
    if(lverb)write(iulog,*)'Constant pft data'
    if(IGO==0)then
      pft_changed_loc=.true.
      iyear=1    
      call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
      call InitPlantMgmnt(NHW,NHE,NVN,NVS)
    elseif(IGO>0)then
      iyear=2
      call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
    endif
   
  elseif(pft_dflag==1)then
    pft_changed_loc=.true.
    if(lverb)write(iulog,*)'Transient pft data',yeari,yearc
    iyear=1
    DO while(.true.)
      call ncd_getvar(pftinfo_nfid,'year',iyear,year)
      if(year==yeari)exit   !when year is found matching the forcing data yeari.
      iyear=iyear+1
    ENDDO
    call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
    call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  else
    pft_changed_loc=.true.
    !transient pft data (flexible inputs)
    iyear=1
    DO while(.true.)
      call ncd_getvar(pftinfo_nfid,'year',iyear,year)
      if(year==yearc)exit
      iyear=iyear+1
    ENDDO
    call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
    call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)  
  endif
  call ncd_pio_closefile(pftinfo_nfid)

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        iHarvestDay_pft(NZ,NY,NX)  = iDayPlantHarvest_pft(NZ,NY,NX)
        iHarvestYear_pft(NZ,NY,NX) = iYearPlantHarvest_pft(NZ,NY,NX)
      ENDDO
    ENDDO
  ENDDO
  if(present(pft_changed))pft_changed=pft_changed_loc
  end subroutine ReadPlantManagementNC
!------------------------------------------------------------------------------------------

  subroutine ReadPlantProperties(nu_plt,NZ,NY,NX,pft_changed)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  logical, intent(in) :: pft_changed
  integer :: N, NB
  real(r8) :: VRNXI,VRNLI

! begin_execution

!
! READ INPUTS FOR EACH PLANT SPECIES
!
  IF(DATAP(NZ,NY,NX).NE.'NO')THEN

    call SetPlantTraits(nu_plt,pft_changed,NZ,NY,NX,VRNLI,VRNXI)
!
!   RE-CALCULATE PLANT INPUTS IN MODEL UNITS
!
!   LeafSWabsorpty_pft,LeafPARabsorpty_pft=leaf SW,PAR absorbtivity
!
    VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)  = 2.5_r8*VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)
    VmaxRubOxyRef_pft(NZ,NY,NX)      = 2.5_r8*VmaxRubOxyRef_pft(NZ,NY,NX)
    VmaxPEPCarboxyRef_pft(NZ,NY,NX)  = 2.5_r8*VmaxPEPCarboxyRef_pft(NZ,NY,NX)
    SpecLeafChlAct_pft(NZ,NY,NX)     = 2.5_r8*SpecLeafChlAct_pft(NZ,NY,NX)
    LeafSWabsorpty_pft(NZ,NY,NX)     = 1.0_r8-RadSWLeafAlbedo_pft(NZ,NY,NX)-RadSWLeafTransmis_pft(NZ,NY,NX)
    LeafPARabsorpty_pft(NZ,NY,NX)    = 1.0_r8-CanopyPARalbedo_pft(NZ,NY,NX)-RadPARLeafTransmis_pft(NZ,NY,NX)
    RadSWLeafAlbedo_pft(NZ,NY,NX)    = RadSWLeafAlbedo_pft(NZ,NY,NX)/LeafSWabsorpty_pft(NZ,NY,NX)
    CanopyPARalbedo_pft(NZ,NY,NX)    = CanopyPARalbedo_pft(NZ,NY,NX)/LeafPARabsorpty_pft(NZ,NY,NX)
    RadSWLeafTransmis_pft(NZ,NY,NX)  = RadSWLeafTransmis_pft(NZ,NY,NX)/LeafSWabsorpty_pft(NZ,NY,NX)
    RadPARLeafTransmis_pft(NZ,NY,NX) = RadPARLeafTransmis_pft(NZ,NY,NX)/LeafPARabsorpty_pft(NZ,NY,NX)
    SineBranchAngle_pft(NZ,NY,NX)    = SIN(BranchAngle_pft(NZ,NY,NX)*RadianPerDegree)
    SinePetioleAngle_pft(NZ,NY,NX)   = SIN(PetioleAngle_pft(NZ,NY,NX)*RadianPerDegree)
    MatureGroup_pft(NZ,NY,NX)        = GROUPX_pft(NZ,NY,NX)

    IF(iPlantTurnoverPattern_pft(NZ,NY,NX).NE.0)THEN
      RefNodeInitRate_pft(NZ,NY,NX)        = RefNodeInitRate_pft(NZ,NY,NX)/MaxNodesPerBranch
      RefLeafAppearRate_pft(NZ,NY,NX)      = RefLeafAppearRate_pft(NZ,NY,NX)/MaxNodesPerBranch
      MatureGroup_pft(NZ,NY,NX)            = MatureGroup_pft(NZ,NY,NX)/MaxNodesPerBranch
      ShootNodeNumAtPlanting_pft(NZ,NY,NX) = ShootNodeNumAtPlanting_pft(NZ,NY,NX)/MaxNodesPerBranch
    ENDIF
    MatureGroup_pft(NZ,NY,NX)=MatureGroup_pft(NZ,NY,NX)-ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    IF(CriticPhotoPeriod_pft(NZ,NY,NX).LT.0.0_r8)THEN
      CriticPhotoPeriod_pft(NZ,NY,NX)=DayLenthMax_col(NY,NX)
    ENDIF
    D5: DO NB=1,NumCanopyLayers
      IF(iPlantPhenolType_pft(NZ,NY,NX).EQ.iphenotyp_evgreen .AND. iPlantPhenolPattern_pft(NZ,NY,NX).NE.iplt_annual)THEN
        HourReq4LeafOut_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNLI+144.0_r8*PlantInitThermoAdaptZone_pft(NZ,NY,NX)*(NB-1))
        HourReq4LeafOff_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNXI+144.0_r8*PlantInitThermoAdaptZone_pft(NZ,NY,NX)*(NB-1))
      ELSE
        HourReq4LeafOut_brch(NB,NZ,NY,NX)=VRNLI
        HourReq4LeafOff_brch(NB,NZ,NY,NX)=VRNXI
      ENDIF
    ENDDO D5
  ENDIF
! WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,CriticPhotoPeriod_pft(NZ,NY,NX)
!1111    FORMAT(A20,2I8,E12.4)
  end subroutine ReadPlantProperties

!------------------------------------------------------------------------------------------
  subroutine ReadPlantTraitsNC(nu_plt,NZ,NY,NX,VRNLI,VRNXI)

  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  real(r8), intent(out) :: VRNLI,VRNXI
  character(len=40) :: pft_lname
  character(len=64):: koppen_climl
  character(len=3) :: koppen_clims
  integer :: loc,N

  loc=get_pft_loc(KoppenClimZone_col(NY,NX),DATAP(NZ,NY,NX)(1:6),pft_lname,koppen_climl,koppen_clims)
  DATAPI(NZ,NY,NX)=loc
  call ncd_getvar(pft_nfid, 'ICTYP', loc, iPlantPhotosynthesisType(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IGTYP', loc, iPlantRootProfile_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ISTYP', loc, iPlantPhenolPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IDTYP', loc, iPlantDevelopPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'INTYP', loc, iPlantNfixType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IWTYP', loc, iPlantPhenolType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IPTYP', loc, iPlantPhotoperiodType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IBTYP', loc, iPlantTurnoverPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IRTYP', loc, iPlantGrainType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'MY', loc, Myco_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ZTYPI', loc, PlantInitThermoAdaptZone_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'VCMX', loc, VmaxSpecRubCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VOMX', loc, VmaxRubOxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VCMX4', loc, VmaxPEPCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO2', loc, XKCO2_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKO2', loc, XKO2_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO24', loc, Km4PEPCarboxy_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RUBP', loc, LeafRuBPConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PEPC', loc, FracLeafProtAsPEPCarboxyl_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ETMX', loc, SpecLeafChlAct_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL', loc, LeafC3ChlorofilConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL4', loc, LeafC4ChlorofilConc_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'FCO2', loc, CanopyCi2CaRatio_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'ALBR', loc, RadSWLeafAlbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ALBP', loc, CanopyPARalbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUR', loc, RadSWLeafTransmis_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUP', loc, RadPARLeafTransmis_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'XRNI', loc, RefNodeInitRate_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XRLA', loc, RefLeafAppearRate_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CTC', loc, TCChill4Seed_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VRNLI', loc,VRNLI)
  call ncd_getvar(pft_nfid, 'VRNXI', loc,VRNXI)
  call ncd_getvar(pft_nfid, 'WDLF', loc,rLen2WidthLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PB', loc,NonstCMinConc2InitBranch_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'GROUPX', loc,GROUPX_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XTLI', loc,ShootNodeNumAtPlanting_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XDL', loc,CriticPhotoPeriod_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XPPD', loc,PhotoPeriodSens_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'SLA1', loc,SLA1_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SSL1', loc,PetoLen2Mass_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SNL1', loc,NodeLenPergC_pft(NZ,NY,NX))


  call ncd_getvar(pft_nfid, 'CLASS', loc,LeafAngleClass_pft(1:NumLeafZenithSectors,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CFI', loc,ClumpFactorInit_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGBR', loc,BranchAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGSH', loc,PetioleAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'STMX', loc,MaxPotentSeedNumber_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SDMX', loc,MaxSeedNumPerSite_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRMX', loc,SeedCMassMax_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRDM', loc,SeedCMass_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GFILL', loc,GrainFillRate25C_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'WTSTDI', loc,StandingDeadInitC_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'RRAD1M', loc,Root1stMaxRadius_pft(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RRAD2M', loc,Root2ndMaxRadius_pft(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PORT', loc,RootPorosity_pft(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PR', loc,MinNonstC2InitRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRR', loc,RootRadialResist_pft(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRA', loc,RootAxialResist_pft(1,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PTSHT', loc,ShutRutNonstElmntConducts_pft(NZ,NY,NX))
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

  call ncd_getvar(pft_nfid, 'OSMO', loc,CanOsmoPsi0pt_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RCS', loc,RCS_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSMX', loc,CuticleResist_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'DMLF', loc,LeafBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSHE', loc,PetioleBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSTK', loc,StalkBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRSV', loc,ReserveBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMHSK', loc,HuskBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMEAR', loc,EarBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMGR', loc,GrainBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRT', loc,RootBiomGrosYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMND', loc,NoduGrowthYield_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CNLF', loc,rNCLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSHE', loc,rNCSheath_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSTK', loc,rNCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRSV', loc,rNCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNHSK', loc,rNCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNEAR', loc,rNCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNGR', loc,rNCGrain_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRT', loc,rNCRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNND', loc,rNCNodule_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CPLF', loc,rPCLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSHE', loc,rPCSheath_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSTK', loc,rPCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRSV', loc,rPCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPHSK', loc,rPCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPEAR', loc,rPCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPGR', loc,rPCGrain_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRT', loc,rPCRootr_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPND', loc,rPCNoduler_pft(NZ,NY,NX))

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
  subroutine SetPlantTraits(nu_plt,pft_changed,NZ,NY,NX,VRNLI,VRNXI)

  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  logical, intent(in) :: pft_changed
  real(r8), intent(out) :: VRNLI,VRNXI

  character(len=40) :: pft_lname
  character(len=64):: koppen_climl
  character(len=3) :: koppen_clims
  integer :: loc,N
  
  loc=get_pft_loc(KoppenClimZone_col(NY,NX),DATAP(NZ,NY,NX)(1:6),pft_lname,koppen_climl,koppen_clims)

  iPlantPhotosynthesisType(NZ,NY,NX)  = iPlantPhotosynthesisType_tab(loc)
  iPlantRootProfile_pft(NZ,NY,NX)     = iPlantRootProfile_tab(loc)
  iPlantPhenolPattern_pft(NZ,NY,NX)   = iPlantPhenolPattern_tab(loc)
  iPlantDevelopPattern_pft(NZ,NY,NX)  = iPlantDevelopPattern_tab(loc)
  iPlantNfixType_pft(NZ,NY,NX)        = iPlantNfixType_tab(loc)
  iPlantPhenolType_pft(NZ,NY,NX)      = iPlantPhenolType_tab(loc)
  iPlantPhotoperiodType_pft(NZ,NY,NX) = iPlantPhotoperiodType_tab(loc)
  iPlantTurnoverPattern_pft(NZ,NY,NX) = iPlantTurnoverPattern_tab(loc)
  iPlantGrainType_pft(NZ,NY,NX)       = iPlantGrainType_tab(loc)
  Myco_pft(NZ,NY,NX)                  = Myco_tab(loc)
  PlantInitThermoAdaptZone_pft(NZ,NY,NX)  = PlantInitThermoAdaptZone_pft_tab(loc)

  VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)         = VmaxRubCarboxyRef_tab(loc)
  VmaxRubOxyRef_pft(NZ,NY,NX)             = VmaxRubOxyRef_tab(loc)
  VmaxPEPCarboxyRef_pft(NZ,NY,NX)         = VmaxPEPCarboxyRef_tab(loc)
  XKCO2_pft(NZ,NY,NX)                     = XKCO2_tab(loc)
  XKO2_pft(NZ,NY,NX)                      = XKO2_tab(loc)
  Km4PEPCarboxy_pft(NZ,NY,NX)             = Km4PEPCarboxy_tab(loc)
  LeafRuBPConc_pft(NZ,NY,NX)              = LeafRuBPConc_tab(loc)
  FracLeafProtAsPEPCarboxyl_pft(NZ,NY,NX) = FracLeafProtAsPEPCarboxyl_tab(loc)
  SpecLeafChlAct_pft(NZ,NY,NX)          = SpecChloryfilAct_tab(loc)
  LeafC3ChlorofilConc_pft(NZ,NY,NX)       = LeafC3ChlorofilConc_tab(loc)
  LeafC4ChlorofilConc_pft(NZ,NY,NX)       = LeafC4ChlorofilConc_tab(loc)
  CanopyCi2CaRatio_pft(NZ,NY,NX)                = CanopyCi2CaRatio_pft_tab(loc)

  RadSWLeafAlbedo_pft(NZ,NY,NX)    = RadSWLeafAlbedo_tab(loc)
  CanopyPARalbedo_pft(NZ,NY,NX)    = CanopyPARalbedo_tab(loc)
  RadSWLeafTransmis_pft(NZ,NY,NX)  = RadSWLeafTransmis_tab(loc)
  RadPARLeafTransmis_pft(NZ,NY,NX) = RadPARLeafTransmis_tab(loc)

  RefNodeInitRate_pft(NZ,NY,NX)          = RefNodeInitRate_tab(loc)
  RefLeafAppearRate_pft(NZ,NY,NX)        = RefLeafAppearRate_tab(loc)
  TCChill4Seed_pft(NZ,NY,NX)             = TCChill4Seed_tab(loc)
  VRNLI                                  = VRNLI_tab(loc)
  VRNXI                                  = VRNXI_tab(loc)
  rLen2WidthLeaf_pft(NZ,NY,NX)           = rLen2WidthLeaf_tab(loc)
  NonstCMinConc2InitBranch_pft(NZ,NY,NX) = NonstCMinConc2InitBranch_tab(loc)

  GROUPX_pft(NZ,NY,NX)                 = GROUPX_tab(loc)
  ShootNodeNumAtPlanting_pft(NZ,NY,NX) = ShootNodeNumAtPlanting_tab(loc)
  CriticPhotoPeriod_pft(NZ,NY,NX)      = CriticPhotoPeriod_tab(loc)
  PhotoPeriodSens_pft(NZ,NY,NX)        = PhotoPeriodSens_tab(loc)

  SLA1_pft(NZ,NY,NX)         = SLA1_tab(loc)
  PetoLen2Mass_pft(NZ,NY,NX) = PetoLen2Mass_tab(loc)
  NodeLenPergC_pft(NZ,NY,NX) = NodeLenPergC_tab(loc)

  LeafAngleClass_pft(1:NumLeafZenithSectors,NZ,NY,NX) = LeafAngleClass_tab(1:NumLeafZenithSectors,loc)
  ClumpFactorInit_pft(NZ,NY,NX)                       = ClumpFactorInit_tab(loc)
  BranchAngle_pft(NZ,NY,NX)                           = BranchAngle_tab(loc)
  PetioleAngle_pft(NZ,NY,NX)                          = PetioleAngle_tab(loc)
  MaxPotentSeedNumber_pft(NZ,NY,NX)                   = MaxPotentSeedNumber_tab(loc)
  MaxSeedNumPerSite_pft(NZ,NY,NX)                     = MaxSeedNumPerSite_tab(loc)
  SeedCMassMax_pft(NZ,NY,NX)                          = SeedCMassMax_tab(loc)
  SeedCMass_pft(NZ,NY,NX)                             = SeedCMass_tab(loc)
  GrainFillRate25C_pft(NZ,NY,NX)                      = GrainFillRate25C_tab(loc)
  StandingDeadInitC_pft(NZ,NY,NX)                     = StandingDeadInitC_tab(loc)

  Root1stMaxRadius_pft(1,NZ,NY,NX)        = Root1stMaxRadius_tab(loc)
  Root2ndMaxRadius_pft(1,NZ,NY,NX)        = Root2ndMaxRadius_tab(loc)
  RootPorosity_pft(1,NZ,NY,NX)            = RootPorosity_tab(loc)
  MinNonstC2InitRoot_pft(NZ,NY,NX)        = MinNonstC2InitRoot_tab(loc)
  RootRadialResist_pft(1,NZ,NY,NX)        = RootRadialResist_tab(loc)
  RootAxialResist_pft(1,NZ,NY,NX)         = RootAxialResist_tab(loc)
  ShutRutNonstElmntConducts_pft(NZ,NY,NX) = ShutRutNonstElmntConducts_tab(loc)
  RootBranchFreq_pft(NZ,NY,NX)            = RootBranchFreq_tab(loc)

  VmaxNH4Root_pft(ipltroot,NZ,NY,NX) = VmaxNH4Root_tab(loc)
  KmNH4Root_pft(ipltroot,NZ,NY,NX)   = KmNH4Root_tab(loc)
  CMinNH4Root_pft(ipltroot,NZ,NY,NX) = CMinNH4Root_tab(loc)

  VmaxNO3Root_pft(ipltroot,NZ,NY,NX) = VmaxNO3Root_tab(loc)
  KmNO3Root_pft(ipltroot,NZ,NY,NX)   = KmNO3Root_tab(loc)
  CminNO3Root_pft(ipltroot,NZ,NY,NX) = CminNO3Root_tab(loc)

  VmaxPO4Root_pft(ipltroot,NZ,NY,NX) = VmaxPO4Root_tab(Loc)
  KmPO4Root_pft(ipltroot,NZ,NY,NX)   = KmPO4Root_tab(loc)
  CMinPO4Root_pft(ipltroot,NZ,NY,NX) = CMinPO4Root_tab(loc)

  CanOsmoPsi0pt_pft(NZ,NY,NX) = CanOsmoPsi0pt_tab(loc)
  RCS_pft(NZ,NY,NX)           = RCS_tab(loc)
  CuticleResist_pft(NZ,NY,NX) = CuticleResist_tab(loc)

  LeafBiomGrowthYld_pft(NZ,NY,NX)    = LeafBiomGrowthYld_tab(loc)
  PetioleBiomGrowthYld_pft(NZ,NY,NX) = PetioleBiomGrowthYld_tab(loc)
  StalkBiomGrowthYld_pft(NZ,NY,NX)   = StalkBiomGrowthYld_tab(loc)
  ReserveBiomGrowthYld_pft(NZ,NY,NX) = ReserveBiomGrowthYld_tab(loc)
  HuskBiomGrowthYld_pft(NZ,NY,NX)    = HuskBiomGrowthYld_tab(loc)
  EarBiomGrowthYld_pft(NZ,NY,NX)     = EarBiomGrowthYld_tab(loc)
  GrainBiomGrowthYld_pft(NZ,NY,NX)   = GrainBiomGrowthYld_tab(loc)
  RootBiomGrosYld_pft(NZ,NY,NX)      = RootBiomGrosYld_tab(loc)
  NoduGrowthYield_pft(NZ,NY,NX)      = NoduGrowthYield_tab(loc)

  rNCLeaf_pft(NZ,NY,NX)    = rNCLeaf_tab(loc)
  rNCSheath_pft(NZ,NY,NX)  = rNCSheath_tab(loc)
  rNCStalk_pft(NZ,NY,NX)   = rNCStalk_tab(loc)
  rNCReserve_pft(NZ,NY,NX) = rNCReserve_tab(loc)
  rNCHusk_pft(NZ,NY,NX)    = rNCHusk_tab(loc)
  rNCEar_pft(NZ,NY,NX)     = rNCEar_tab(loc)
  rNCGrain_pft(NZ,NY,NX)   = rNCGrain_tab(loc)
  rNCRoot_pft(NZ,NY,NX)    = rNCRoot_tab(loc)
  rNCNodule_pft(NZ,NY,NX)  = rNCNodule_tab(loc)

  rPCLeaf_pft(NZ,NY,NX)    = rPCLeaf_tab(loc)
  rPCSheath_pft(NZ,NY,NX)  = rPCSheath_tab(loc)
  rPCStalk_pft(NZ,NY,NX)   = rPCStalk_tab(loc)
  rPCReserve_pft(NZ,NY,NX) = rPCReserve_tab(loc)
  rPCHusk_pft(NZ,NY,NX)    = rPCHusk_tab(loc)
  rPCEar_pft(NZ,NY,NX)     = rPCEar_tab(loc)
  rPCGrain_pft(NZ,NY,NX)   = rPCGrain_tab(loc)
  rPCRootr_pft(NZ,NY,NX)   = rPCRootr_tab(loc)
  rPCNoduler_pft(NZ,NY,NX) = rPCNoduler_tab(loc)

  if(disp_planttrait.and.pft_changed)then
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

  end subroutine SetPlantTraits

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
  character(len=60) :: strval

!   iPlantPhotosynthesisType=photosynthesis type:3=C3,4=C4
!   iPlantRootProfile_pft=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!   iPlantDevelopPattern_pft=growth habit:0=determinate,1=indetermimate
!   iPlantNfixType_pft=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   iPlantPhenolType_pft=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   iPlantPhotoperiodType_pft=photoperiod type:0=day neutral,1=short day,2=long day
!   iPlantTurnoverPattern_pft=turnover:if iPlantRootProfile_pft=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if iPlantRootProfile_pft=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   iPlantGrainType_pft=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   PlantInitThermoAdaptZone_pft=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical
  write(nu_plt,*)('=',j=1,100)
  write(nu_plt,*)'PLANT traits for FUNCTIONAL TYPE (NZ,NY,NX)=',NZ,NY,NX,trim(DATAP(NZ,NY,NX)(1:6))
  call writefixsl(nu_plt,'Plant name ',pft_lname,40)
  strval=koppen_clims//','//koppen_climl
  call writefixsl(nu_plt,'Koppen climate info',strval,40)

  select CASE (iPlantPhotosynthesisType(NZ,NY,NX))
  case (3)
    strval='C3'
  case (4)
    strval='C4'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Photosynthesis pathway',strval,40)

  select case(iPlantRootProfile_pft(NZ,NY,NX))
  case (0)
    strval='Shallow root profile, like bryophytes'
  case (1)
    strval='Intermediate root profile, like herbs'
  case (2)
    strval='Deep root profile, like trees'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Root profile pattern',strval,40)

  select case (iPlantPhenolPattern_pft(NZ,NY,NX))
  case (0)
    strval='Annual'
  case (1)
    strval='Perennial'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Life cycle',strval,40)

  select case (iPlantDevelopPattern_pft(NZ,NY,NX))
  case (0)
    strval='Determinate'
  case (1)
    strval='Indetermimate'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Growth pattern IDTYP',strval,40)

  select case (iPlantNfixType_pft(NZ,NY,NX))
! 1,2, 3, e.g. legumes
  case (in2fixtyp_root_fast)
    strval='Rapid root N-fixation'
  case (in2fixtyp_root_medium)
    strval='Intermediate root N-fixation'
  case (in2fixtyp_root_slow)
    strval='Slow root N-fixation'
!4,5,6, e.g. cyanobacteria
  case (in2fixtyp_canopy_fast)
    strval='Rapid canopy N-fixation'
  case (in2fixtyp_canopy_medium)
    strval='Intermediate canopy N-fixation'
  case (in2fixtyp_canopy_slow)
    strval='Slow canopy N-fixation'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'N-fixation symbiosis INTYP',strval,40)

  select case(iPlantPhenolType_pft(NZ,NY,NX))
  case (iphenotyp_evgreen)
    write(strval,'(A,X,I3)')'Evergreen',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_coldecid)
    write(strval,'(A,X,I3)')'Cold deciduous',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_drouhtdecidu)
    write(strval,'(A,X,I3)')'Drought deciduous',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_coldroutdecid)
    write(strval,'(A,X,I3)')'Cold&drought-tolerant deciduous',iPlantPhenolType_pft(NZ,NY,NX)
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Phenology type IWTYP',strval,40)

  select case(iPlantPhotoperiodType_pft(NZ,NY,NX))
  case (iphotop_neutral)
    strval='Neutral day'
  case (iphotop_short)
    strval='Short day'
  case (iphotop_long)
    strval='Long day'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Photoperiod IPTYP',strval,40)

  if(is_plant_treelike(iPlantRootProfile_pft(NZ,NY,NX)))then
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(strval,'(A,I2)')'Rapid, like deciduous tree ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (2)
      write(strval,'(A,I2)')'Very slow, like evergreen needleleaf tree ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (3)
      write(strval,'(A,I2)')'Slow, like evergreen broadleaf ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (4)
      write(strval,'(A,I2)')'Like semi-deciduous tree ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (5)
      write(strval,'(A,I2)')'Like semi-evergreen tree ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case default
      write(strval,'(A,I2)')'Undefined tree-like pattern ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    end select
  else
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(strval,'(A,I2)')'Rapid all aboveground biome, herbaceous ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case default
      write(strval,'(A,I2)')'Undefined herbaceous pattern ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    end select
  endif
  call writefixsl(nu_plt,'Biome turnover pattern IBTYP',strval,40)

  select case(iPlantGrainType_pft(NZ,NY,NX))
  case (igraintyp_abvgrnd)
    strval='Aboveground 0'
  case (igraintyp_blwgrnd)
    strval='Belowground 1'
  case default
    strval='Not defined'
  end select
  call writefixsl(nu_plt,'Storage organ IRTYP',strval,40)

  select case(Myco_pft(NZ,NY,NX))
  case (1)
    strval='No'
  case (2)
    strval='YES'
  case default
    strval='Wrong option'
  end select
  call writefixsl(nu_plt,'Mycorrhizal association MY',strval,40)

  if(PlantInitThermoAdaptZone_pft(NZ,NY,NX)<0 .or. PlantInitThermoAdaptZone_pft(NZ,NY,NX)>ithermozone_tropical)then
    write(strval,'(A,X,F6.2)')'Not defined',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif(PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_arcboreal)then
    write(strval,'(A,X,F6.2)')'Arctic, boreal',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_cooltempr)then
    write(strval,'(A,X,F6.2)')'Cool temperate',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_warmtempr)then
    write(strval,'(A,X,F6.2)')'Warm temperate',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_subtropic)then
    write(strval,'(A,X,F6.2)')'Subtropical',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_tropical)then
    write(strval,'(A,X,F6.2)')'Tropical',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  endif
  
  call writefixsl(nu_plt,'Thermal adaptation zone ZTYPI',strval,40)
  end subroutine pft_display

!------------------------------------------------------------------------------------------
  subroutine photosyns_trait_disp(nu_plt,NZ,NY,NX)
  !
  !Description
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j
  character(len=81) :: line
  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'PHOTOSYNTHETIC PROPERTIES'
  call writefixl(nu_plt,'Specific rubisco carboxylase [umol C g-1 s-1] VCMX',VmaxSpecRubCarboxyRef_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Specific rubisco oxygenase [umol O2 g-1 s-1] VOMX',VmaxRubOxyRef_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Specific PEP carboxylase activity [umol g-1 s-1] VCMX4',VmaxPEPCarboxyRef_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Km for VmaxSpecRubCarboxyRef_pft [uM] XKCO2',XKCO2_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Km for VmaxRubOxyRef_pft [uM] XKO2',XKO2_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'KM for VmaxPEPCarboxyRef_pft [uM] XKCO24',Km4PEPCarboxy_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Fraction of leaf protein in rubisco [g rub/(gC protein)] RUBP',LeafRuBPConc_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Fraction of leaf protein in PEP carboxylase [g PEP/(gC protein)] PEPC',&
    FracLeafProtAsPEPCarboxyl_pft(NZ,NY,NX),90)
  call writefixl(nu_plt,'Specific chlorophyll activity [umol e- gC-1 s-1] ETMX',SpecLeafChlAct_pft(NZ,NY,NX),90)
  if(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic3_photo)then
    call writefixl(nu_plt,'Fraction of leaf protein as chlorophyll in mesophyll (C3) [g Chl /(gC protein)] CHL',&
      LeafC3ChlorofilConc_pft(NZ,NY,NX),90)
  elseif(iPlantPhotosynthesisType(NZ,NY,NX).eq.ic4_photo)then
    call writefixl(nu_plt,'Fraction of leaf protein as chlorophyll in bundle sheath(C4) [g Chl /(gC protein)] CHL',&
      LeafC3ChlorofilConc_pft(NZ,NY,NX),90)
    call writefixl(nu_plt,'Fraction of leaf protein in mesophyll chlorophyll(C4) [g Chl /(gC protein)] CHL4',&
      LeafC4ChlorofilConc_pft(NZ,NY,NX),90)
  endif
  call writefixl(nu_plt,'Intercellular-to-atmospheric CO2 concentration ratio [-]FCO2',CanopyCi2CaRatio_pft(NZ,NY,NX),90)
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
  call writefixl(nu_plt,'Rate of node initiation at 25oC [h-1] XRNI',RefNodeInitRate_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Rate of leaf appearance at 25oC [h-1] XRLA',RefLeafAppearRate_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Chilling temperature for CO2 fixation, seed loss [oC] CTC',TCChill4Seed_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Hour requirement for spring leafout [h] VRNLI',VRNLI,84)
  call writefixl(nu_plt,'Hour requirement for autumn leafoff [h] VRNXI',VRNXI,84)
  call writefixl(nu_plt,'Leaf length:width ratio WDLF [-]',rLen2WidthLeaf_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Nonstructural C concentration needed for branching PB [gC]',NonstCMinConc2InitBranch_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Maturity group, node number required for floral initiation [-] GROUPX',GROUPX_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Node number at planting [-] XTLI',ShootNodeNumAtPlanting_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Critical photoperiod <= maximum daylength [h] XDL',CriticPhotoPeriod_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Photoperiod sensitivity [node h-1] XPPD',PhotoPeriodSens_pft(NZ,NY,NX),84)


  end subroutine Phenology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine morphology_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: N,j

  write(nu_plt,*)('-',j=1,100)  
  write(nu_plt,*)'MORPHOLOGICAL PROPERTIES'
  call writefixl(nu_plt,'growth in leaf area vs mass (m2 g-1) SLA1',SLA1_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'growth in petiole length vs mass (m g-1) SSL1',PetoLen2Mass_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'growth in internode length vs mass (m g-1) SNL1',NodeLenPergC_pft(NZ,NY,NX),84)
  call writeafixl(nu_plt,'fraction of leaf area in 0-22.5,45,67.5,90o inclination classes CLASS',&
    LeafAngleClass_pft(1:NumLeafZenithSectors,NZ,NY,NX),84)
  call writefixl(nu_plt,'initial clumping factor CFI',ClumpFactorInit_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'stem angle from horizontal ANGBR',BranchAngle_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'petiole angle from horizontal ANGSH',PetioleAngle_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'maximum potential seed mumber from '// &
    'pre-anthesis stalk growth STMX',MaxPotentSeedNumber_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'maximum seed number per STMX (none) SDMX',MaxSeedNumPerSite_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'maximum seed size per SDMX (gC) GRMX',SeedCMassMax_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'seed size at planting (gC) GRDM',SeedCMass_pft(NZ,NY,NX),84)    !could be greater than SeedCMassMax, accouting for seedling
  call writefixl(nu_plt,'grain filling rate at 25 oC (g seed-1 h-1) GFILL',GrainFillRate25C_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'mass of dead standing biomass at planting (gC m-2) WTSTDI',StandingDeadInitC_pft(NZ,NY,NX),84)
  end subroutine morphology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ROOT CHARACTERISTICS'
  call writefixl(nu_plt,'Radius of fine roots [m] RRAD1M',Root1stMaxRadius_pft(1,NZ,NY,NX),84)
  call writefixl(nu_plt,'Radius of root hairs [m] RRAD2M',Root2ndMaxRadius_pft(1,NZ,NY,NX),84)
  call writefixl(nu_plt,'Primary/fine root porosity [m3 m-3] PORT',RootPorosity_pft(1,NZ,NY,NX),84)
  call writefixl(nu_plt,'Nonstructural C concentration needed for root'// &
    ' branching (gC/gC) PR',MinNonstC2InitRoot_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Radial resistivity per m2 root surface area for water uptake [MPa h m-1] RSRR',RootRadialResist_pft(1,NZ,NY,NX),84)
  call writefixl(nu_plt,'Axial resistivity per m root length for water uptake [MPa h m-4] RSRA',RootAxialResist_pft(1,NZ,NY,NX),84)
  call writefixl(nu_plt,'Rate constant for equilibrating shoot-root '// &
    'nonstructural C concn PTSHT',ShutRutNonstElmntConducts_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'Branching frequency from 1st to 2nd roots [m-1] RTFQ',RootBranchFreq_pft(NZ,NY,NX),84)
  end subroutine Root_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ROOT UPTAKE PARAMETERS'
  call writefixl(nu_plt,'NH4 max uptake [g m-2 h-1] UPMXZH',VMaxNH4Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'NH4 uptake Km [uM] UPKMZH',KmNH4Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'NH4 uptake min concn [uM] UPMNZH',CMinNH4Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'NO3 max uptake [g m-2 h-1] UPMXZO',VMaxNO3Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'NO3 uptake Km [uM] UPKMZO',KmNO3Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'NO3 uptake min concn [uM] UPMNZO',CMinNO3Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'H2PO4 or H1PO4 max uptake [g m-2 h-1] UPMXPO',VMaxPO4Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'H2PO4 or H1PO4 uptake Km [uM] UPKMPO',KmPO4Root_pft(ipltroot,NZ,NY,NX),84)
  call writefixl(nu_plt,'H2PO4 or H1PO4 uptake min conc [uM] UPMNPO',CMinPO4Root_pft(ipltroot,NZ,NY,NX),84)
  end subroutine Root_nutrient_trait_disp


!------------------------------------------------------------------------------------------
  subroutine plant_water_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)  
  write(nu_plt,*)'WATER RELATIONS'
  call writefixl(nu_plt,'leaf osmotic potential at zero leaf water '// &
     'potential [MPa] OSMO',CanOsmoPsi0pt_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'shape parameter for stomatal resistance vs '// &
     'leaf turgor potential [MPa-1] RCS',RCS_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'H2O cuticular resistance [s m-1] RSMX',CuticleResist_pft(NZ,NY,NX),84)
  end subroutine plant_water_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ORGAN GROWTH YIELDS'
  call writefixl(nu_plt,'leaf dry matter C production vs '// &
    'nonstructural C consumption (g g-1) DMLF',LeafBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'petiole dry matter C production vs '// &
    'nonstructural C consumption (g g-1) DMSHE',PetioleBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'stalk dry matter C production vs '// &
    'nonstructural C consumption (g g-1) DMSTK',StalkBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'stalk reserve C production vs '// &
    'nonstructural C consumption (g g-1) DMRSV',ReserveBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'husk dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1) DMHSK',HuskBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'ear dry matter C production vs '// &
    'nonstructural Cconsumption (g g-1) DMEAR',EarBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'grain C production vs nonstructural C'// &
    ' consumption (g g-1) DMGR',GrainBiomGrowthYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'root dry matter C production vs nonstructural C'// &
    ' consumption (g g-1) DMRT',RootBiomGrosYld_pft(NZ,NY,NX),101)
  call writefixl(nu_plt,'nodule bacteria in root, canopy dry matter '// &
    'C production vs nonstructural C consumption (g g-1) DMND' &
    ,NoduGrowthYield_pft(NZ,NY,NX),101)
  end subroutine plant_biomyield_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'ORGAN N AND P CONCENTRATIONS'
  call writefixl(nu_plt,'NC ratio in plant leaves (gN/gC) CNLF',rNCLeaf_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant petiole (gN/gC) CNSHE',rNCSheath_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant stalk (gN/gC) CNSTK',rNCStalk_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant stalk reserve (gN/gC) CNRSV',rNCReserve_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant husk (gN/gC) CNHSK',rNCHusk_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant ear (gN/gC) CNEAR',rNCEar_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant grain (gN/gC) CNGR',rNCGrain_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant root (gN/gC) CNRT',rNCRoot_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'NC ratio in plant nodule (gN/gC) CNND',rNCNodule_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant leaves (gP/gC) CPLF',rPCLeaf_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant petiole (gP/gC) CPSHE',rPCSheath_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant stalk (gP/gC) CPSTK',rPCStalk_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant stalk reserve (gP/gC) CPRSV',rPCReserve_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant husk (gP/gC) CPHSK',rPCHusk_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant ear (gP/gC) CPEAR',rPCEar_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant grain (gP/gC) CPGR',rPCGrain_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant root (gP/gC) CPRT',rPCRootr_pft(NZ,NY,NX),84)
  call writefixl(nu_plt,'PC ratio in plant nodule (gP/gC) CPND',rPCNoduler_pft(NZ,NY,NX),84)
  end subroutine plant_biomstoich_trait_disp

!------------------------------------------------------------------------------------------

  subroutine plant_optic_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j

  write(nu_plt,*)('-',j=1,100)
  write(nu_plt,*)'OPTICAL PROPERTIES'
  call writefixl(nu_plt,'leaf SW albedo ALBR',RadSWLeafAlbedo_pft(NZ,NY,NX),50)
  call writefixl(nu_plt,'leaf PAR albedo ALBP',CanopyPARalbedo_pft(NZ,NY,NX),50)
  call writefixl(nu_plt,'leaf SW transmission TAUR',RadSWLeafTransmis_pft(NZ,NY,NX),50)
  call writefixl(nu_plt,'leaf PAR transmission TAUP',RadPARLeafTransmis_pft(NZ,NY,NX),50)
  end subroutine plant_optic_trait_disp


!------------------------------------------------------------------------------------------

  subroutine ReadPlantInfoNC(yeari,NHW,NHE,NVN,NVS)
  !
  !Description
  use netcdf
  use ncdio_pio
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: yeari    !forcing year information
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
  integer :: ret,kk
  character(len=1) :: cc
  character(len=10) :: pft_gtype(JP)

  IF(is_first_year)THEN
    D9999: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        DO NZ=1,JP
          iDayPlanting_pft(NZ,NY,NX)      = -1E+06
          iYearPlanting_pft(NZ,NY,NX)     = -1E+06
          iDayPlantHarvest_pft(NZ,NY,NX)  = 1E+06
          iYearPlantHarvest_pft(NZ,NY,NX) = 1E+06
        enddo
      enddo
    ENDDO D9999

  ENDIF

  if (len_trim(pft_mgmt_in)==0)then
   !no plant management info to read
    D5995: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NP_col(NY,NX)=0
        DO NZ=1,NP0_col(NY,NX)
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
          NP_col(NY,NX)=0
          DO NZ=1,NP0_col(NY,NX)
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
      if(first_topou)ntopou=1    
      DO NTOPO=1,ntopou
        call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
        call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
        call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
        call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)

        NN=1
        if(any((/NH1, NV1,NH2,NV2/)<=0) .or. any((/NH1, NV1,NH2,NV2/)>ntopou))exit
        D4995: DO NX=NH1,NH2
          D4990: DO NY=NV1,NV2
            NP0_col(NY,NX)=NS
            D4985: DO NZ=1,NS
              LSG_pft(NZ,NY,NX)=NN
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
                flag_active_pft(NZ,NY,NX)=.true.
                !modify pft name 
                IF(KoppenClimZone_col(NY,NX).GT.0)THEN
                  WRITE(CLIMATE,'(I2)')KoppenClimZone_col(NY,NX)
                  !the type of pft is specified by genra+Koppen climate zone
                  DATAX_pft(NZ)=trim(pft_gtype(NZ)(1:4)//CLIMATE)                  
                  !consider cases when koppen climate zone is zero  
                elseif(KoppenClimZone_col(NY,NX)==0)then
                  DATAX_pft(NZ)=''
                  DO kk=1,10
                    cc=pft_gtype(NZ)(kk:kk)
                    if(ichar(cc)==0)exit
                    DATAX_pft(NZ)(kk:kk)=cc
                  ENDDO
                ENDIF
              ENDDO D4965
              if(first_pft)then
                flag_active_pft(2:NS,NY,NX)=.false.
                NP0_col(NY,NX)=1
                NS=1
              endif  
            ENDDO D4970
          ENDDO D4975
        ENDIF

        D8995: DO NX=NH1,NH2
          D8990: DO NY=NV1,NV2
            NP_col(NY,NX)=NS
            !DATAP(NZ,NY,NX) and DATAM(NZ,NY,NX) are to be read in readqmod.F90
            D100: DO NZ=1,NP_col(NY,NX)
              DATAP(NZ,NY,NX)=DATAX_pft(NZ)
            ENDDO D100

            D101: DO NZ=NP_col(NY,NX)+1,JP
              DATAP(NZ,NY,NX)='NO'
            ENDDO D101
          ENDDO D8990
        ENDDO D8995
      ENDDO
    ENDIF 
    call ncd_pio_closefile(pftinfo_nfid)
  endif
  end subroutine ReadPlantInfoNC

!------------------------------------------------------------------------------------------

  subroutine writefixl(nu_plt,desc,value,width)
  implicit none
  integer,intent(in)  :: nu_plt
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  real(r8),intent(in) :: value

  character(len=width) :: line
  line=desc
  line(width:width)=':'
  write(nu_plt,*)line,value
  end subroutine writefixl
!------------------------------------------------------------------------------------------

  subroutine writefixsl(nu_plt,desc,value,width)
  implicit none
  integer,intent(in)  :: nu_plt
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  character(len=*),intent(in) :: value

  character(len=width) :: line
  line=desc
  line(width:width)=':'
  write(nu_plt,*)line,' ',value
  end subroutine writefixsl

!------------------------------------------------------------------------------------------

  subroutine writeafixl(nu_plt,desc,values,width)
  implicit none
  integer,intent(in)  :: nu_plt
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  real(r8),dimension(:),intent(in) :: values

  character(len=width) :: line
  line=desc
  line(width:width)=':'
  write(nu_plt,*)line,values
  end subroutine writeafixl
!------------------------------------------------------------------------------------------

  Subroutine ReadPlantTraitTable()

  implicit none

  call ncd_pio_openfile(pft_nfid, pft_file_in, ncd_nowrite)
  call ncd_getvar(pft_nfid, 'pfts', pftss_tab)
  call ncd_getvar(pft_nfid,'pfts_long',pft_long_tab)
  call ncd_getvar(pft_nfid,'pfts_short',pft_short_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_no',koppen_clim_ncode_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_short',koppen_clim_short_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_long',koppen_clim_long_tab)
  call ncd_getvar(pft_nfid, 'ICTYP', iPlantPhotosynthesisType_tab)
  call ncd_getvar(pft_nfid, 'IGTYP', iPlantRootProfile_tab)
  call ncd_getvar(pft_nfid, 'ISTYP', iPlantPhenolPattern_tab)
  call ncd_getvar(pft_nfid, 'IDTYP', iPlantDevelopPattern_tab)
  call ncd_getvar(pft_nfid, 'INTYP', iPlantNfixType_tab)
  call ncd_getvar(pft_nfid, 'IWTYP', iPlantPhenolType_tab)
  call ncd_getvar(pft_nfid, 'IPTYP', iPlantPhotoperiodType_tab)
  call ncd_getvar(pft_nfid, 'IBTYP', iPlantTurnoverPattern_tab)
  call ncd_getvar(pft_nfid, 'IRTYP', iPlantGrainType_tab)
  call ncd_getvar(pft_nfid, 'MY',  Myco_tab)
  call ncd_getvar(pft_nfid, 'ZTYPI', PlantInitThermoAdaptZone_pft_tab)
  call ncd_getvar(pft_nfid, 'VCMX',  VmaxRubCarboxyRef_tab)
  call ncd_getvar(pft_nfid, 'VOMX',  VmaxRubOxyRef_tab)
  call ncd_getvar(pft_nfid, 'VCMX4', VmaxPEPCarboxyRef_tab)
  call ncd_getvar(pft_nfid, 'XKCO2', XKCO2_tab)
  call ncd_getvar(pft_nfid, 'XKO2',  XKO2_tab)
  call ncd_getvar(pft_nfid, 'XKCO24',Km4PEPCarboxy_tab)
  call ncd_getvar(pft_nfid, 'RUBP', LeafRuBPConc_tab)
  call ncd_getvar(pft_nfid, 'PEPC', FracLeafProtAsPEPCarboxyl_tab)
  call ncd_getvar(pft_nfid, 'ETMX', SpecChloryfilAct_tab)
  call ncd_getvar(pft_nfid, 'CHL',  LeafC3ChlorofilConc_tab)
  call ncd_getvar(pft_nfid, 'CHL4', LeafC4ChlorofilConc_tab)
  call ncd_getvar(pft_nfid, 'FCO2', CanopyCi2CaRatio_pft_tab)
  call ncd_getvar(pft_nfid, 'ALBR', RadSWLeafAlbedo_tab)
  call ncd_getvar(pft_nfid, 'ALBP', CanopyPARalbedo_tab)
  call ncd_getvar(pft_nfid, 'TAUR', RadSWLeafTransmis_tab)
  call ncd_getvar(pft_nfid, 'TAUP', RadPARLeafTransmis_tab)
  call ncd_getvar(pft_nfid, 'XRNI', RefNodeInitRate_tab)
  call ncd_getvar(pft_nfid, 'XRLA', RefLeafAppearRate_tab)
  call ncd_getvar(pft_nfid, 'CTC',  TCChill4Seed_tab)
  call ncd_getvar(pft_nfid, 'VRNLI', VRNLI_tab)
  call ncd_getvar(pft_nfid, 'VRNXI', VRNXI_tab)
  call ncd_getvar(pft_nfid, 'WDLF', rLen2WidthLeaf_tab)
  call ncd_getvar(pft_nfid, 'PB', NonstCMinConc2InitBranch_tab)
  call ncd_getvar(pft_nfid, 'GROUPX', GROUPX_tab)
  call ncd_getvar(pft_nfid, 'XTLI', ShootNodeNumAtPlanting_tab)
  call ncd_getvar(pft_nfid, 'XDL', CriticPhotoPeriod_tab)
  call ncd_getvar(pft_nfid, 'XPPD',PhotoPeriodSens_tab)
  call ncd_getvar(pft_nfid, 'SLA1', SLA1_tab)
  call ncd_getvar(pft_nfid, 'SSL1', PetoLen2Mass_tab)
  call ncd_getvar(pft_nfid, 'SNL1', NodeLenPergC_tab)
  call ncd_getvar(pft_nfid, 'CLASS', LeafAngleClass_tab)
  call ncd_getvar(pft_nfid, 'CFI', ClumpFactorInit_tab)
  call ncd_getvar(pft_nfid, 'ANGBR', BranchAngle_tab)
  call ncd_getvar(pft_nfid, 'ANGSH', PetioleAngle_tab)
  call ncd_getvar(pft_nfid, 'STMX', MaxPotentSeedNumber_tab)
  call ncd_getvar(pft_nfid, 'SDMX', MaxSeedNumPerSite_tab)
  call ncd_getvar(pft_nfid, 'GRMX', SeedCMassMax_tab)
  call ncd_getvar(pft_nfid, 'GRDM', SeedCMass_tab)
  call ncd_getvar(pft_nfid, 'GFILL',GrainFillRate25C_tab)
  call ncd_getvar(pft_nfid, 'WTSTDI',StandingDeadInitC_tab)
  call ncd_getvar(pft_nfid, 'RRAD1M',Root1stMaxRadius_tab)
  call ncd_getvar(pft_nfid, 'RRAD2M',Root2ndMaxRadius_tab)
  call ncd_getvar(pft_nfid, 'PORT', RootPorosity_tab)
  call ncd_getvar(pft_nfid, 'PR', MinNonstC2InitRoot_tab)
  call ncd_getvar(pft_nfid, 'RSRR', RootRadialResist_tab)
  call ncd_getvar(pft_nfid, 'RSRA', RootAxialResist_tab)
  call ncd_getvar(pft_nfid, 'PTSHT',ShutRutNonstElmntConducts_tab)
  call ncd_getvar(pft_nfid, 'RTFQ', RootBranchFreq_tab)
  call ncd_getvar(pft_nfid, 'UPMXZH',VmaxNH4Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMZH',KmNH4Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNZH',CMinNH4Root_tab)
  call ncd_getvar(pft_nfid, 'UPMXZO', VmaxNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMZO', KmNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNZO', CminNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPMXPO', VmaxPO4Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMPO', KmPO4Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNPO', CMinPO4Root_tab)
  call ncd_getvar(pft_nfid, 'OSMO', CanOsmoPsi0pt_tab)
  call ncd_getvar(pft_nfid, 'RCS', RCS_tab)
  call ncd_getvar(pft_nfid, 'RSMX',CuticleResist_tab)
  call ncd_getvar(pft_nfid, 'DMLF', LeafBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMSHE',PetioleBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMSTK',StalkBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMRSV',ReserveBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMHSK',HuskBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMEAR',EarBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMGR', GrainBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMRT', RootBiomGrosYld_tab)
  call ncd_getvar(pft_nfid, 'DMND', NoduGrowthYield_tab)
  call ncd_getvar(pft_nfid, 'CNLF', rNCLeaf_tab)
  call ncd_getvar(pft_nfid, 'CNSHE', rNCSheath_tab)
  call ncd_getvar(pft_nfid, 'CNSTK', rNCStalk_tab)
  call ncd_getvar(pft_nfid, 'CNRSV', rNCReserve_tab)
  call ncd_getvar(pft_nfid, 'CNHSK', rNCHusk_tab)
  call ncd_getvar(pft_nfid, 'CNEAR', rNCEar_tab)
  call ncd_getvar(pft_nfid, 'CNGR', rNCGrain_tab)
  call ncd_getvar(pft_nfid, 'CNRT', rNCRoot_tab)
  call ncd_getvar(pft_nfid, 'CNND', rNCNodule_tab)
  call ncd_getvar(pft_nfid, 'CPLF', rPCLeaf_tab)
  call ncd_getvar(pft_nfid, 'CPSHE',rPCSheath_tab)
  call ncd_getvar(pft_nfid, 'CPSTK', rPCStalk_tab)
  call ncd_getvar(pft_nfid, 'CPRSV', rPCReserve_tab)
  call ncd_getvar(pft_nfid, 'CPHSK', rPCHusk_tab)
  call ncd_getvar(pft_nfid, 'CPEAR', rPCEar_tab)
  call ncd_getvar(pft_nfid, 'CPGR', rPCGrain_tab)
  call ncd_getvar(pft_nfid, 'CPRT', rPCRootr_tab)
  call ncd_getvar(pft_nfid, 'CPND', rPCNoduler_tab)
  call ncd_pio_closefile(pft_nfid)
  end Subroutine ReadPlantTraitTable

!------------------------------------------------------------------------------------------

  function get_pft_loc(koppen_def,pft_name,pft_lname,koppen_climl,koppen_clims)result(loc)
  !
  !!DESCRIPTION
  ! return the id of pft to be read
  use PlantTraitTableMod
  implicit none
  integer, intent(in) :: koppen_def            !koppen climate numerical code
  character(len=*), intent(in) :: pft_name     !pft short name
  character(len=40),intent(out):: pft_lname    !returning pft long name
  character(len=64),intent(out):: koppen_climl !koppen climate Description
  character(len=3), intent(out):: koppen_clims !koppen climate short code
  integer :: loc,loc1,len,k

  len=len_trim(pft_name)

  do 
    if(ichar(pft_name(len:len))==0 .or. pft_name(len:len)==' ')then
      len=len-1
    else
      exit
    endif
  enddo

  loc=1
  DO
    if(pftss_tab(loc)(1:len)==pft_name(1:len))exit
    loc=loc+1
  enddo

  loc1=1
  do 
    if(pft_name(1:4)==pft_short_tab(loc1))exit
    loc1=loc1+1
    if(loc1 > size(pft_short_tab))then
      exit      
    endif
  enddo
  if(loc1 <=  size(pft_short_tab))then
    pft_lname=pft_long_tab(loc1)
  else
    pft_lname=pft_name
  endif
  
  if(koppen_def==0)then
    koppen_climl='None'  
    koppen_clims='None'
    return
  endif
  loc1=1  
  do
    if(koppen_clim_ncode_tab(loc1)==pft_name(5:6))exit
    loc1=loc1+1
  enddo
  koppen_climl=koppen_clim_long_tab(loc1)
  koppen_clims=koppen_clim_short_tab(loc1)
  end function get_pft_loc
end module PlantInfoMod
