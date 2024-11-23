module RestartMod
!
! DESCRIPTION
! code to write restart/check point files
! code to read/write restart files
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use data_const_mod    , only : spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL  
  use EcoSIMConfig      , only : jcplx=> jcplxc, NumMicbFunGrupsPerCmplx=> NumMicbFunGrupsPerCmplx,nlbiomcp=>NumLiveMicrbCompts
  use EcoSIMConfig      , only : ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc, cold_run
  use EcoSIMConfig      , only : inst_suffix,ref_date,start_date, ctitle, finidat
  use EcoSIMConfig      , only : case_name,hostname,version,source,username
  use EcoSiMParDataMod  , only : micpar,pltpar
  use TracerIDMod       , only : trc_confs
  use EcoSIMCtrlMod     , only : etimer,do_budgets,plant_model
  use ecosim_log_mod   , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun,destroy
  use HistFileMod       , only : hist_restart_ncd  
  use EcoSIMCtrlDataType  
  use restUtilMod  
  use ncdio_pio
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use CanopyRadDataType
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SnowDataType
  use SurfLitterDataType
  use SurfSoilDataType
  use CanopyDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use AqueChemDatatype
  use SedimentDataType
  use GridDataType
  use fileUtil
implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  character(len=8),public :: rest_opt='year'
  integer         ,public :: rest_frq=-999999999
  integer,  parameter :: uninit_int = -999999999
  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), parameter :: namet  = 'topounit'     ! name of topographic units
  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches

  public :: restFile
  public :: get_restart_date
  contains

!------------------------------------------------------------------------------------------
  subroutine restFile(flag)
  implicit none
  character(len=*), intent(in) :: flag
  character(len=256) :: fnamer       ! name of netcdf restart file
  character(len=256) :: filer                   ! restart file name
  character(len=256) :: path       ! name of netcdf restart file  
  character(len=18) :: rdate
  integer :: yr,mon,day,tod

  if (flag=='read')then
    call restFile_getfile(fnamer,path)
    call restFile_read( bounds, fnamer)    
  else if(flag=='write')then  
    call etimer%get_curr_date(yr,mon,day,tod)
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i6.6)') yr,mon,day,tod
    filer = restFile_filename(rdate=rdate)
    call restFile_write(bounds, filer,rdate=rdate)
  endif
  end subroutine restFile

!------------------------------------------------------------------------------------------
  SUBROUTINE restartnc_plant(ncid,flag)  
!
!     THIS SUBROUTINE WRITES OUT ALL PLANT VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id

  integer :: NHW,NHE,NVN,NVS
  integer :: K,L,NX,NY,NZ,N,NB,NR,M
  real(r8), target :: datrc_1d(bounds%begc:bounds%endc)
  integer, target  :: datic_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrp_1d(bounds%begp:bounds%endp)
  real(r8), target :: datrp_2d(bounds%begp:bounds%endp,1:MaxNumBranches+1)
  integer , target :: datip_1d(bounds%begp:bounds%endp)
  integer , target :: datip_2d(bounds%begp:bounds%endp,1:MaxNumBranches+1)
  real(r8), target, allocatable :: datrp_3d(:,:,:)
  real(r8), target, allocatable :: datrp_4d(:,:,:,:)  
  real(r8), target, allocatable :: datrp_5d(:,:,:,:,:)  
  integer , target, allocatable :: datip_3d(:,:,:)
  integer, pointer :: dat1pr(:),dat2pr(:,:),dat3pr(:,:,:)
  real(r8),pointer :: datpr1(:),datpr2(:,:),datpr3(:,:,:)
  real(r8),pointer :: datpr4(:,:,:,:),datpr5(:,:,:,:,:)
  integer :: ncols, npfts
  integer :: ic,ip,sz3,sz4,sz5,sz2

! execution begins here
  NHW = bounds%NHW;NVN = bounds%NVN
  NHE = bounds%NHE;NVS = bounds%NVS
  ncols = bounds%ncols
  npfts = bounds%npfts
  sz2=max(MaxNumBranches+1,JZ+1,NumOfCanopyLayers)
  sz3=max(jsken,NumPlantChemElms,JZ+1,MaxNodesPerBranch+1)
  sz4=max(MaxNodesPerBranch+1,NumOfCanopyLayers)
  sz5=max(NumOfLeafZenithSectors,NumPlantChemElms)

  allocate(datrp_3d(bounds%begp:bounds%endp,sz3,sz2))
  allocate(datrp_4d(bounds%begp:bounds%endp,sz4,sz3,sz2))  
  allocate(datrp_5d(bounds%begp:bounds%endp,sz5,sz4,sz3,sz2))  
  allocate(datip_3d(bounds%begp:bounds%endp,sz3,MaxNumBranches+1))

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NP', dim1name='column',&
       long_name='number of plant species', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NP,datic_1d)
  else 
    !print*,'NP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NP,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NP', dim1name='column',&
       long_name='number of plant species', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='LeafStalkArea_col', dim1name='column',&
       long_name='leaf stalk area', units='m2 d-2', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LeafStalkArea_col,datrc_1d)
  else     
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LeafStalkArea_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='LeafStalkArea_col', dim1name='column',&
       long_name='leaf stalk area', units='m2 d-2', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RoughHeight', dim1name='column',&
       long_name='Roughness height', units='m', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,RoughHeight_col,datrc_1d)
  else     
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RoughHeight_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RoughHeight', dim1name='column',&
       long_name='Roughness height', units='m', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
  endif


  if(flag=='read')then
    dat1pr=>datic_1d
    call restartvar(ncid, flag, varname='NumActivePlants', dim1name='column',&
       long_name='number of active PFT', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NumActivePlants,datic_1d)

  else 
    !print*,'NumActivePlants'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NumActivePlants,datic_1d)
    dat1pr=>datic_1d
    call restartvar(ncid, flag, varname='NumActivePlants', dim1name='column',&
       long_name='number of active PFT', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='DATAPI', dim1name='pft',&
       long_name='ID for pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,DATAPI,datip_1d)
  else 
    !print*,'DATAPI'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,DATAPI,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='DATAPI', dim1name='pft',&
       long_name='ID for pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='IsPlantActive_pft', dim1name='pft',&
       long_name='flag for living pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,IsPlantActive_pft,datip_1d)

  else 
    !print*,'IsPlantActive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,IsPlantActive_pft,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='IsPlantActive_pft', dim1name='pft',&
       long_name='flag for living pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iPlantState_pft', dim1name='pft',&
       long_name='flag for species death, 0:live, 1:dead', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantState_pft,datip_1d)
  else 
    !print*,'iPlantState_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantState_pft,datip_1d)  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iPlantState_pft', dim1name='pft',&
       long_name='flag for species death, 0:live, 1:dead', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='NumCogrothNode_pft', dim1name='pft',&
       long_name='Number of growing plant nodes', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumCogrothNode_pft,datip_1d)
  else 
    !print*,'iPlantState_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NumCogrothNode_pft,datip_1d)  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='NumCogrothNode_pft', dim1name='pft',&
       long_name='Number of growing plant nodes', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:jsken)
    call restartvar(ncid, flag, varname='StandDeadKCompElms_pft', dim1name='pft',dim2name='nkinecmp',&
       dim3name='elmnts',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadKCompElms_pft,datrp_3d)

  else
    !print*,'StandDeadKCompElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadKCompElms_pft,datrp_3d)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:jsken)    
    call restartvar(ncid, flag, varname='StandDeadKCompElms_pft', dim1name='pft',dim2name='elmnts',&
       dim3name='nkinecmp',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)

  endif

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlanting_pft', dim1name='pft',&
       long_name='year of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'iYearPlanting_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlanting_pft', dim1name='pft',&
       long_name='year of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iDayPlanting_pft', dim1name='pft',&
       long_name='day of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'iDayPlanting_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iDayPlanting_pft', dim1name='pft',&
       long_name='day of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlantHarvest_pft', dim1name='pft',&
       long_name='year of harvest', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)

  else
    !print*,'iYearPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlantHarvest_pft,datip_1d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iYearPlantHarvest_pft', dim1name='pft',&
       long_name='year of harvest', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    

  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iDayPlantHarvest_pft', dim1name='pft',&
     long_name='day of harvest', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else 
    !print*,'iDayPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlantHarvest_pft,datip_1d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iDayPlantHarvest_pft', dim1name='pft',&
     long_name='day of harvest', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif 

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NGTopRootLayer_pft', dim1name='pft',&
     long_name='top soil layer that has root', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NGTopRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'NGTopRootLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NGTopRootLayer_pft,datip_1d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NGTopRootLayer_pft', dim1name='pft',&
     long_name='top soil layer that has root', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iPlantShootState_pft', dim1name='pft',&
     long_name='flag to detect canopy death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantShootState_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)

  else
    !print*,'iPlantShootState_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantShootState_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iPlantShootState_pft', dim1name='pft',&
     long_name='flag to detect canopy death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)      
  endif

  if(flag=='read')then    
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iPlantRootState_pft', dim1name='pft',&
     long_name='flag to detect root system death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantRootState_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'iPlantRootState_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantRootState_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iPlantRootState_pft', dim1name='pft',&
     long_name='flag to detect root system death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        

  endif

  if(flag=='read')then   
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumOfBranches_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumOfBranches_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else 
    !print*,'NumOfBranches_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NumOfBranches_pft,datip_1d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumOfBranches_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)            
  endif 

  if(flag=='read')then     
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='BranchNumber_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else 
    !print*,'BranchNumber_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='BranchNumber_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)            
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MainBranchNum_pft', dim1name='pft',&
     long_name='number of main branch', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MainBranchNum_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else 
    !print*,'MainBranchNum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,MainBranchNum_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MainBranchNum_pft', dim1name='pft',&
     long_name='number of main branch', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)     

  endif     

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='doInitPlant_pft', dim1name='pft',&
     long_name='PFT initialization:0=no,1=yes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitPlant_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'doInitPlant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitPlant_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='doInitPlant_pft', dim1name='pft',&
     long_name='PFT initialization:0=no,1=yes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))         
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MaxSoiL4Root_pft', dim1name='pft',&
     long_name='maximum root layer', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MaxSoiL4Root_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,MaxSoiL4Root_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MaxSoiL4Root_pft', dim1name='pft',&
     long_name='maximum root layer', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))         
  endif  

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumRootAxes_pft', dim1name='pft',&
     long_name='root primary axis number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumRootAxes_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'NRT'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NumRootAxes_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumRootAxes_pft', dim1name='pft',&
     long_name='root primary axis number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NIXBotRootLayer_pft', dim1name='pft',&
     long_name='Bottom soil layer number that has root', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else 
    !print*,'NIXBotRootLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NIXBotRootLayer_pft', dim1name='pft',&
     long_name='Bottom soil layer number that has root', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)     

  endif 

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MY', dim1name='pft',&
     long_name='mycorrhizal type:0-no,1-yes', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
     fill_value=ispval,flag_values=(/0,1/))       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MY,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'MY'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,MY,datip_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MY', dim1name='pft',&
     long_name='mycorrhizal type:0-no,1-yes', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
     fill_value=ispval,flag_values=(/0,1/))       
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)    
    call restartvar(ncid, flag, varname='NIXBotRootLayer_rpft', dim1name='pft',dim2name='rootaxs',&
     long_name='maximum soil layer number for root axes', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, &
     fill_value=ispval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_rpft,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'NIXBotRootLayer_rpft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_rpft,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)  
    call restartvar(ncid, flag, varname='NIXBotRootLayer_rpft', dim1name='pft',dim2name='rootaxs',&
     long_name='maximum soil layer number for root axes', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)     
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)  
    call restartvar(ncid, flag, varname='PlantExudElm_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantExudElm_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PlantExudElm_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantExudElm_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)  
    call restartvar(ncid, flag, varname='PlantExudElm_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='LitrfalStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element LitrFall', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LitrfalStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'LitrfalStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LitrfalStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='LitrfalStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element LitrFall', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantN2Fix_CumYr_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantN2Fix_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PlantN2Fix_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantN2Fix_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantN2Fix_CumYr_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_N_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootUptk_N_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootUptk_N_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_N_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_P_CumYr_pft', dim1name='pft',&
     long_name='total plant P uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootUptk_P_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootUptk_P_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_P_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossResp_pft', dim1name='pft',&
     long_name='total plant respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossResp_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'GrossResp_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossResp_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossResp_pft', dim1name='pft',&
     long_name='total plant respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CNRTS_pft', dim1name='pft',&
     long_name='root N:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CNRTS_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CNRTS_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CNRTS_pft', dim1name='pft',&
     long_name='root N:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CPRTS_pft', dim1name='pft',&
     long_name='root P:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CPRTS_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CPRTS_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CPRTS_pft', dim1name='pft',&
     long_name='root P:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Uptk_pft', dim1name='pft',&
     long_name='threshold for plant nutrient uptake', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4Uptk_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4Uptk_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Uptk_pft', dim1name='pft',&
     long_name='threshold for plant nutrient uptake', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantO2Stress_pft', dim1name='pft',&
     long_name='plant O2 stress', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantO2Stress_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  

  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantO2Stress_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantO2Stress_pft', dim1name='pft',&
     long_name='plant O2 stress', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4LeafVar_pft', dim1name='pft',&
     long_name='threshold for plant leaf dynamics', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4LeafVar_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4LeafVar_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4LeafVar_pft', dim1name='pft',&
     long_name='threshold for plant leaf dynamics', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Groth_pft', dim1name='pft',&
     long_name='threshold for plant growth', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4Groth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ZERO4Groth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Groth_pft', dim1name='pft',&
     long_name='threshold for plant growth', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif


  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2Canopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,FracPARads2Canopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,FracPARads2Canopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2Canopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedTempSens_pft', dim1name='pft',&
     long_name='seed temperature sensitivity', units='1/oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedTempSens_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedTempSens_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedTempSens_pft', dim1name='pft',&
     long_name='seed temperature sensitivity', units='1/oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ETCanopy_CumYr_pft', dim1name='pft',&
     long_name='total transpiration (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ETCanopy_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'ETCanopy_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ETCanopy_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ETCanopy_CumYr_pft', dim1name='pft',&
     long_name='total transpiration  (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossCO2Fix_pft', dim1name='pft',&
     long_name='total gross CO2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossCO2Fix_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'GrossCO2Fix_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossCO2Fix_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossCO2Fix_pft', dim1name='pft',&
     long_name='total gross CO2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCelciusCanopy_pft', dim1name='pft',&
     long_name='canopy temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TCelciusCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TCelciusCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TCelciusCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCelciusCanopy_pft', dim1name='pft',&
     long_name='canopy temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='H2OCuticleResist_pft', dim1name='pft',&
     long_name='Cutitle resistance for H2O', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,H2OCuticleResist_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,H2OCuticleResist_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='H2OCuticleResist_pft', dim1name='pft',&
     long_name='Cutitle resistance for H2O', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootFracRemobilizableBiom', dim1name='pft',&
     long_name='fraction of remobilizable nonstructural biomass in root', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootFracRemobilizableBiom,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootFracRemobilizableBiom,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootFracRemobilizableBiom', dim1name='pft',&
     long_name='fraction of remobilizable nonstructural biomass in root', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKC', dim1name='pft',&
     long_name='canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TKC,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TKC'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TKC,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKC', dim1name='pft',&
     long_name='canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TCGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TCG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TCGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TKGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TKG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TKGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           

  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTCanopyGroth_pft', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    call cppft(flag,NHW,NHE,NVN,NVS,NP,fTCanopyGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'fTCanopyGroth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,fTCanopyGroth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTCanopyGroth_pft', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStalkC_pft', dim1name='pft',&
     long_name='canopy active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStalkC_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyStalkC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStalkC_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStalkC_pft', dim1name='pft',&
     long_name='canopy active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               

  endif

  if(flag=='read')then    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyWater_pft', dim1name='pft',&
     long_name='plant canopy water content in dry matter', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyWater_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyWater_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyWater_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyWater_pft', dim1name='pft',&
     long_name='canopy water content in dry matter', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopy_pft', dim1name='pft',&
     long_name='canopy total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PSICanP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopy_pft', dim1name='pft',&
     long_name='canopy total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyOsmo_pft', dim1name='pft',&
     long_name='plant canopy osmotic water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyOsmo_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PSICanopyOsmo_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyOsmo_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyOsmo_pft', dim1name='pft',&
     long_name='plant canopy osmotic water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyTurg_pft', dim1name='pft',&
     long_name='plant canopy turgor water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyTurg_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PSICanopyTurg_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyTurg_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyTurg_pft', dim1name='pft',&
     long_name='plant canopy turgor water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2CuticleResist_pft', dim1name='pft',&
     long_name='CO2 cuticle resistance', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2CuticleResist_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PSICanopyTurg_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2CuticleResist_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2CuticleResist_pft', dim1name='pft',&
     long_name='CO2 cuticle resistance', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootBiomCPerPlant_pft', dim1name='pft',&
     long_name='root C per plant', units='g g-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootBiomCPerPlant_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  

  else
    !print*,'RootBiomCPerPlant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootBiomCPerPlant_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootBiomCPerPlant_pft', dim1name='pft',&
     long_name='root C per plant', units='g g-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafArea_pft', dim1name='pft',&
     long_name='plant leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafArea_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyLeafArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafArea_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafArea_pft', dim1name='pft',&
     long_name='plant leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStemArea_pft', dim1name='pft',&
     long_name='plant stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyStemArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStemArea_pft', dim1name='pft',&
     long_name='plant stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopulation_pft', dim1name='pft',&
     long_name='plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantPopulation_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantPopulation_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopulation_pft', dim1name='pft',&
     long_name='plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracGroth2Node_pft', dim1name='pft',&
     long_name='parameter for allocation of growth to nodes', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,FracGroth2Node_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,FracGroth2Node_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracGroth2Node_pft', dim1name='pft',&
     long_name='parameter for allocation of growth to nodes', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmnt_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'EcoHavstElmnt_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmnt_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant shoot elements', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'EcoHavstElmnt_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant shoot elements', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootMycoExudElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root exudation', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMycoExudElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMycoExudElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootMycoExudElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root exudation', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmntCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
     call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmntCum_pft,datrp_2d,NumActivePlants=NumActivePlants,&
       IsPlantActive_pft=IsPlantActive_pft)  
  else 
    !print*,'EcoHavstElmntCum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmntCum_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmntCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyRespC_CumYr_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyRespC_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyRespC_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyRespC_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyRespC_CumYr_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NetCumElmntFlx2Plant_pft', dim1name='pft',dim2name='elmnts',&
     long_name='effect of canopy element status on seed set', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NetCumElmntFlx2Plant_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'NetCumElmntFlx2Plant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NetCumElmntFlx2Plant_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NetCumElmntFlx2Plant_pft', dim1name='pft',dim2name='elmnts',&
     long_name='effect of canopy element status on seed set', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3Emis_CumYr_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3Emis_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'NH3Emis_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3Emis_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3Emis_CumYr_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SurfLitrfalStrutElms_CumYr_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SurfLitrfalStrutElms_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LitrfalStrutElms_CumYr_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LitrfalStrutElms_CumYr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI_pft',dim1name='pft',&
     long_name='initial plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PPI_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PPI'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PPI_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI_pft',dim1name='pft',&
     long_name='initial plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX_pft', dim1name='pft',&
     long_name='plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PPX_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PPX'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PPX_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX_pft', dim1name='pft',&
     long_name='plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='standing dead element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'StandDeadStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='standing dead element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyHeight_pft', dim1name='pft',&
     long_name='canopy height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyHeight_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyHeight'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyHeight_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyHeight_pft', dim1name='pft',&
     long_name='canopy height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='WatByPCanopy_pft', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            
    call cppft(flag,NHW,NHE,NVN,NVS,NP,WatByPCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'WatByPCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,WatByPCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='WatByPCanopy_pft', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  
  
  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor_pft', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ClumpFactor_pft,datrp_1d,NumActivePlants=NumActivePlants &
      ,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'ClumpFactor'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ClumpFactor_pft,datrp_1d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor_pft', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
  else
    !print*,'CO2ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CH4ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
  else  !define and write
    !print*,'CH4ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CH4ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ENGYX_pft', dim1name='pft',&
     long_name='canopy heat storage from previous time step', units='MJ d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ENGYX_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ENGYX_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ENGYX_pft', dim1name='pft',&
     long_name='canopy heat storage from previous time step', units='MJ d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,O2ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'O2ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,O2ByFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3byFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'NH3byFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3byFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,N2ObyFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'N2ObyFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,N2ObyFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PO4byFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'PO4byFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PO4byFire_CumYr_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  
  
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='dReproNodeNumNormByMatG_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant repro node normalized by maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,dReproNodeNumNormByMatG_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'GROUP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,dReproNodeNumNormByMatG_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='dReproNodeNumNormByMatG_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant repro node normalized by maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='MatureGroup_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MatureGroup_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'GROUP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,MatureGroup_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='MatureGroup_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNum_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootNodeNum_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'ShootNodeNum_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootNodeNum_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNum_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNum2InitFloral_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNum2InitFloral_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'NodeNum2InitFloral_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNum2InitFloral_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNum2InitFloral_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumberAtAnthesis_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumberAtAnthesis_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'NodeNumberAtAnthesis_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumberAtAnthesis_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumberAtAnthesis_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NumOfLeaves_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumOfLeaves_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
  else
    !print*,'NumOfLeaves_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NumOfLeaves_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NumOfLeaves_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafNumberAtFloralInit_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafNumberAtFloralInit_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
  else
    !print*,'VSTGX'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafNumberAtFloralInit_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafNumberAtFloralInit_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
  else
    !print*,'NodeNumNormByMatgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ReprodNodeNumNormByMatrgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'ReprodNodeNumNormByMatrgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ReprodNodeNumNormByMatrgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotalNodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TotalNodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'TotalNodeNumNormByMatgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TotalNodeNumNormByMatgrp_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotalNodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotReproNodeNumNormByMatrgrp_brch', &
     dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TotReproNodeNumNormByMatrgrp_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'TotReproNodeNumNormByMatrgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TotReproNodeNumNormByMatrgrp_brch,&
      datrp_2d,NumActivePlants=NumActivePlants, &
      IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotReproNodeNumNormByMatrgrp_brch', &
      dim1name='pft',dim2name='nbranches',&
      long_name='normalized node number during reproductive growth stages', units='none', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4Leafout_brch', dim1name='pft',dim2name='nbranches',&
     long_name='heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4Leafout_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4Leafout_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4Leafout_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4Leafout_brch', dim1name='pft',dim2name='nbranches',&
     long_name='heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LeafOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LeafOff_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4LeafOff_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LeafOff_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LeafOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LenthenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LenthenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4LenthenPhotoPeriod_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LenthenPhotoPeriod_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LenthenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4ShortenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4ShortenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4ShortenPhotoPeriod_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4ShortenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4ShortenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours2LeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours2LeafOut_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours2LeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours2LeafOut_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours2LeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HourFailGrainFill_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect physiological maturity from  grain fill', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HourFailGrainFill_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'HourFailGrainFill_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HourFailGrainFill_brch,&
      datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HourFailGrainFill_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect physiological maturity from  grain fill', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HoursDoingRemob_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursDoingRemob_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'FLGZ'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursDoingRemob_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HoursDoingRemob_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitLeafOut_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doInitLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitLeafOut_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)            
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)

  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeafOut_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doPlantLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeafOut_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doPlantLeaveOff_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Prep4Literfall_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Prep4Literfall_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Prep4Literfall_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)            
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4LiterfalAftMature_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantBranchState_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantBranchState_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'iPlantBranchState_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantBranchState_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantBranchState_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'BranchNumber_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KLeafNumber_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'KLEAF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KLeafNumber_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'KVSTG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'KLowestGroLeafNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat3pr => datip_3d(1:npfts,1:NumGrowthStages,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantCalendar_brch', dim1name='pft',dim2name='ngrstages', &
     dim3name='nbranches',long_name='plant growth stage', units='none', &
     interpinic_flag='skip', data=dat3pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantCalendar_brch,datip_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'iPlantCalendar'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantCalendar_brch,datip_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    dat3pr => datip_3d(1:npfts,1:NumGrowthStages,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantCalendar_brch', dim1name='pft',dim2name='ngrstages', &
     dim3name='nbranches',long_name='plant growth stage', units='none', &
     interpinic_flag='skip', data=dat3pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyNonstElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNonstElms_pft,&
      datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif 

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nodule nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyNodulNonstElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nodule nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='VHeatCapCanP_pft', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,VHeatCapCanP_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'VHeatCapCanP_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,VHeatCapCanP_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='VHeatCapCanP_pft', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DeltaTKC_pft', dim1name='pft',&
     long_name='change in canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,DeltaTKC_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'DeltaTKC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,DeltaTKC_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DeltaTKC_pft', dim1name='pft',&
     long_name='change in canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy shoot element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'ShootStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy shoot element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafShethC_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafShethC_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyLeafShethC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafShethC_pft,datrp_1d, &
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafShethC_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='LeafStalkArea_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStalkArea_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafStalkArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStalkArea_pft,datrp_1d, &
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='LeafStalkArea_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='MatureGroup_pft', dim1name='pft',&
     long_name='Plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MatureGroup_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafStalkArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,MatureGroup_pft,datrp_1d, &
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='MatureGroup_pft', dim1name='pft',&
     long_name='Plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanPDailyMin', dim1name='pft',&
     long_name='Daily minimum canopy water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanPDailyMin,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafStalkArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanPDailyMin,datrp_1d, &
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanPDailyMin', dim1name='pft',&
     long_name='Daily minimum canopy water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='AvgCanopyBiomC2Graze_pft', dim1name='pft',&
     long_name='Mean canopy C for grazing', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,AvgCanopyBiomC2Graze_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafStalkArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,AvgCanopyBiomC2Graze_pft,datrp_1d, &
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='AvgCanopyBiomC2Graze_pft', dim1name='pft',&
     long_name='Mean canopy C for grazing', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'RootElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LeafStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'LeafStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LeafStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetoleStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'PetoleStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetoleStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'StalkStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkRsrvElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'StalkRsrvElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkRsrvElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='HuskStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'HuskStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='HuskStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EarStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'EarStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EarStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='GrainStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'GrainStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='GrainStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NodulStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='root total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodulStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'NodulStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NodulStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NodulStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='root total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SeasonalNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant stored nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeasonalNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants &
      ,IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'NonStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeasonalNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SeasonalNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant stored nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HypoctoHeight_pft', dim1name='pft',&
     long_name='Hypocotyledon height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HypoctoHeight_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'HypoctoHeight_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HypoctoHeight_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HypoctoHeight_pft', dim1name='pft',&
     long_name='Hypocotyledon height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedDepth_pft', dim1name='pft',&
     long_name='seeding depth', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedDepth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedDepth_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedDepth_pft', dim1name='pft',&
     long_name='seeding depth', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedMeanLen_pft', dim1name='pft',&
     long_name='mean seed length', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedMeanLen_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedMeanLen_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedMeanLen_pft', dim1name='pft',&
     long_name='mean seed length', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedVolumeMean_pft', dim1name='pft',&
     long_name='mean seed volume', units='m3', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedVolumeMean_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedVolumeMean_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedVolumeMean_pft', dim1name='pft',&
     long_name='mean seed volume', units='m3', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedAreaMean_pft', dim1name='pft',&
     long_name='mean seed area', units='m2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedAreaMean_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedAreaMean_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedAreaMean_pft', dim1name='pft',&
     long_name='mean seed area', units='m2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOut_pft', dim1name='pft',&
     long_name='leaf out temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TC4LeafOut_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TC4LeafOut_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOut_pft', dim1name='pft',&
     long_name='leaf out temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOff_pft', dim1name='pft',&
     long_name='leaf off temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TC4LeafOff_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TC4LeafOff_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOff_pft', dim1name='pft',&
     long_name='leaf off temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HoursTooLowPsiCan_pft', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY)', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursTooLowPsiCan_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
  else
    !print*,'HoursTooLowPsiCan_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursTooLowPsiCan_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HoursTooLowPsiCan_pft', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY)', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ChillHours_pft', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ChillHours_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
  else
    !print*,'CHILL'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ChillHours_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ChillHours_pft', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
  else
    !print*,'RootStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)         
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkBiomassC_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'StalkBiomassC_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkBiomassC_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafPetolBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafPetolBiomassC_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafPetolBiomassC_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafPetolBiomassC_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafPetolBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyNonstElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule nonstructural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyNodulNonstElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule nonstructural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch shoot element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'ShootStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch shoot element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyNodulStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PetoleStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'StalkStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkRsrvElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'StalkRsrvElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkRsrvElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HuskStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'HuskStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HuskStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EarStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'EarStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EarStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'GrainStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',&
      dim2name='nbranches',long_name='branch seed sites', units=' d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PotentialSeedSites_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PotentialSeedSites_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PotentialSeedSites_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',&
      dim2name='nbranches',long_name='branch seed sites', units=' d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SeedNumSet_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch grain number', units='# d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedNumSet_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'SeedNumSet_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedNumSet_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SeedNumSet_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch grain number', units='# d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainSeedBiomCMean_brch', dim1name='pft',dim2name='nbranches',&
     long_name='maximum grain C during grain fill', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainSeedBiomCMean_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'GrainSeedBiomCMean_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainSeedBiomCMean_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainSeedBiomCMean_brch', dim1name='pft',dim2name='nbranches',&
     long_name='maximum grain C during grain fill', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaLive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaLive_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafAreaLive_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaLive_brch,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaLive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafChemElmRemob_brch', dim1name='pft',dim2name='elmnts', &
     dim3name='nbranches',long_name='branch leaf structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafChemElmRemob_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafChemElmRemob_brch', dim1name='pft',dim2name='elmnts', &
     dim3name='nbranches',long_name='branch leaf structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaDying_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaDying_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafAreaDying_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaDying_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaDying_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanPBranchHeight', dim1name='pft',dim2name='nbranches',&
     long_name='branch height', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanPBranchHeight,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanPBranchHeight'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanPBranchHeight,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanPBranchHeight', dim1name='pft',dim2name='nbranches',&
     long_name='branch height', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='element translocated from leaf during senescence', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafElmntRemobFlx_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='element translocated from leaf during senescence', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleChemElmRemob_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch sheath structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PetioleChemElmRemob_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleChemElmRemob_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch sheath structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleChemElmRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='element translocated from sheath during senescence', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PetioleChemElmRemobFlx_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleChemElmRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='element translocated from sheath during senescence', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SenecStalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch stalk structural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SenecStalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'SenecStalkStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SenecStalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SenecStalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch stalk structural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafAreaZ_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyLeafAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafAreaZ_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafCLyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafCLyr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyLeafAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafCLyr_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafCLyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemAreaZ_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyStemAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemAreaZ_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL4_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C4 content in C4 photosynthesi', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL4_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CPOOL4_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL4_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL4_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C4 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL3_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C3 content in C3 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL3_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CPOOL3_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL3_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL3_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C3 content in C3 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassCO2BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassCO2BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CMassCO2BundleSheath_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassCO2BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassCO2BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassHCO3BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassHCO3BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CMassHCO3BundleSheath_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassHCO3BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassHCO3BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'ARLF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafElmntNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinCNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafProteinCNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinCNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleLensNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleLensNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PetoleLensNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleLensNode_brch,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleLensNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif 

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'WGSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)      
  endif  
 
  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleProteinCNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'WSSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleProteinCNode_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LiveInterNodeHight_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LiveInterNodeHight_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LiveInterNodeHight_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LiveInterNodeHight_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LiveInterNodeHight_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='living internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightDying_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightDying_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'InternodeHeightDying_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightDying_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightDying_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeStrutElms_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'InternodeStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeStrutElms_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumOfCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafArea_lpft', dim1name='pft',dim2name='levcan',&
     dim3name='nodes1',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafArea_lpft,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 

  else
    !print*,'CanopyLeafArea_lpft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafArea_lpft,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumOfCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafArea_lpft', dim1name='pft',dim2name='levcan',&
     dim3name='nodes1',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumOfCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmsByLayerNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='levcan',dim4name='nodes1',dim5name='nbranches',long_name='layer leaf element',&
     units='g d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmsByLayerNode_brch,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafElmsByLayerNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmsByLayerNode_brch,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumOfCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmsByLayerNode_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='levcan',dim4name='nodes1',dim5name='nbranches',long_name='layer leaf element',&
      units='g d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr5 => datrp_5d(1:npfts,1:NumOfLeafZenithSectors,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='leaf surface area',&
      units='m2 d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaZsec_brch,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafAreaZsec_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaZsec_brch,datrp_5d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
    datpr5 => datrp_5d(1:npfts,1:NumOfLeafZenithSectors,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='leaf surface area',&
      units='m2 d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:NumOfCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyStalkArea_lbrch', dim1name='pft',dim2name='levcan',&
     dim3name='nbranches',long_name='plant canopy layer branch stem layer area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStalkArea_lbrch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyStalkArea_lbrch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStalkArea_lbrch,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumOfCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyStalkArea_lbrch', dim1name='pft',dim2name='levcan',&
     dim3name='nbranches',long_name='plant canopy layer branch stem layer area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumOfLeafZenithSectors,1:NumOfCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StemAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nbranches',long_name='stem surface area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StemAreaZsec_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'StemAreaZsec_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StemAreaZsec_brch,datrp_4d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumOfLeafZenithSectors,1:NumOfCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StemAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nbranches',long_name='stem surface area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
  endif  
  if(plant_model)then
  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  
  endif
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootVolPerMassC_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root volume per C mass', units='m3 g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootVolPerMassC_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootVolPerMassC_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootVolPerMassC_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root volume per C mass', units='m3 g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootPoreTortu4Gas', dim1name='pft',dim2name='rootyps',&
     long_name='root pore tortuosity for gas transp', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootPoreTortu4Gas,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootPoreTortu4Gas,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootPoreTortu4Gas', dim1name='pft',dim2name='rootyps',&
     long_name='root pore tortuosity for gas transp', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXSecArea_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXSecArea_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stXSecArea_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stXSecArea_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  
  
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stMaxRadius1_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stMaxRadius1_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndMaxRadius1_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndMaxRadius1_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootRaidus_rpft', dim1name='pft',dim2name='rootyps',&
     long_name='root internal radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootRaidus_rpft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootRaidus_rpft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootRaidus_rpft', dim1name='pft',dim2name='rootyps',&
     long_name='root internal radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndSpecLen_pft', dim1name='pft',dim2name='rootyps',&
     long_name='specific root length secondary axes', units='m g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndSpecLen_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndSpecLen_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndSpecLen_pft', dim1name='pft',dim2name='rootyps',&
     long_name='specific root length secondary axes', units='m g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(plant_model)then
    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootOSMO_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root osmotic water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootOSMO_vr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'PSIRootOSMO_vr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootOSMO_vr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootOSMO_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root osmotic water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  
  
    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootTurg_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root turgor water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootTurg_vr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'PSIRootTurg_vr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootTurg_vr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootTurg_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root turgor water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stXNumL_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root/myco number primary axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stXNumL_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root1stXNumL_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stXNumL_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stXNumL_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root/myco number primary axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndXNum_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer 2nd order root number axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root2ndXNum_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndXNum_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer 2nd order root number axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)        
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootLenPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenDensPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root length density', units='m m-3', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenDensPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootLenDensPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenDensPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenDensPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root length density', units='m m-3', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootPoreVol_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume air', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootPoreVol_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootPoreVol_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootPoreVol_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootPoreVol_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume air', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootVH2O_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume water', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootVH2O_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RTVLW'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootVH2O_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootVH2O_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume water', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stRadius_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root1stRadius_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stRadius_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndRadius_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root2ndRadius_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndRadius_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootAreaPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer area per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootAreaPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootAreaPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootAreaPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootAreaPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer area per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 

    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndAveLen_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndAveLen_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root2ndAveLen_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndAveLen_pvr,datrp_3d,&
        NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndAveLen_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootRespPotent_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration unconstrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootRespPotent_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootRespPotent_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootRespPotent_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootRespPotent_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration unconstrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2Autor_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootCO2Autor_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootCO2Autor_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootCO2Autor_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2Autor_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2EmisPot_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootCO2EmisPot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootCO2EmisPot_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootCO2EmisPot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2EmisPot_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then 
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcg_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
      call cppft(flag,NHW,NHE,NVN,NVS,NP,trcg_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'trcg_rootml_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,trcg_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcg_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
    endif  

    if(flag=='read')then 
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcs_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
      call cppft(flag,NHW,NHE,NVN,NVS,NP,trcs_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'trcs_rootml_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,trcs_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcs_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoActiveBiomC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer structural C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMycoActiveBiomC_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootMycoActiveBiomC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMycoActiveBiomC_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoActiveBiomC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer structural C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP, PopuRootMycoC_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,' PopuRootMycoC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP, PopuRootMycoC_pvr,datrp_3d,&
        NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer protein C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootProteinC_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootProteinC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootProteinC_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer protein C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif 

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootO2Dmnd4Resp_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootO2Dmnd4Resp_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootO2Dmnd4Resp_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootO2Dmnd4Resp_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootO2Dmnd4Resp_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNH4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNH4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNH4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNO3DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNO3DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNO3DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH2PO4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootH2PO4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH2PO4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH1PO4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootH1PO4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH1PO4DmndSoil_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNH4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNH4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNH4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNO3DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNO3DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNO3DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH2PO4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootH2PO4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH2PO4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH1PO4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootH1PO4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootH1PO4DmndBand_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RAutoRootO2Limter_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RAutoRootO2Limter_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'WFR'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RAutoRootO2Limter_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RAutoRootO2Limter_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP, RootMycoNonstElms_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,' RootMycoNonstElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP, RootMycoNonstElms_rpvr,datrp_4d,&
        NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root protein concentration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootProteinConc_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'PSIRoot'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootProteinConc_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root protein concentration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNonstructElmConc_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element concentration', units='none', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP, RootNonstructElmConc_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP, RootNonstructElmConc_rpvr,datrp_4d,&
        NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNonstructElmConc_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element concentration', units='none', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:pltpar%jroots,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root1stDepz_pft', dim1name='pft',dim2name='rootyps',&
      dim3name='levcan',long_name='primary root layer depth', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stDepz_pft,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root1stDepz_pft'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stDepz_pft,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts, 1:pltpar%jroots,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root1stDepz_pft', dim1name='pft',dim2name='rootyps',&
      dim3name='levcan',long_name='primary root layer depth', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stElm_raxs', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stElm_raxs,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootMyco1stElm_raxs'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stElm_raxs,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stElm_raxs', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes number', units='d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RTN2'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes number', units='d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root1stLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer length primary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stLen_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root1stLen_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stLen_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root1stLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer length primary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='levcan',long_name='root layer primary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootMyco1stStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='levcan',long_name='root layer primary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='levcan',long_name='root layer secondary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootMyco2ndStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='levcan',long_name='root layer secondary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root2ndLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer length secondary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndLen_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'Root2ndLen_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndLen_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
      call restartvar(ncid, flag, varname='Root2ndLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='levcan',long_name='root layer length secondary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulNonstElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulNonstElms_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNodulNonstElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulNonstElms_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)  
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulNonstElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulStrutElms_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft) 
    else
      !print*,'RootNodulStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulStrutElms_rpvr,datrp_3d,NumActivePlants=NumActivePlants,&
        IsPlantActive_pft=IsPlantActive_pft)   
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  
  endif

  call destroy(datrp_3d)
  call destroy(datrp_4d)  
  call destroy(datrp_5d)  
  call destroy(datip_3d)

  END SUBROUTINE restartnc_plant

!------------------------------------------------------------------------------------------

  subroutine restartnc(ncid,flag)
  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id

  call restartnc_plant(ncid,flag)

  call restartnc_soil(ncid,flag)

  end subroutine restartnc
!------------------------------------------------------------------------------------------

  SUBROUTINE restartnc_soil(ncid,flag)
!
!     THIS SUBROUTINE WRITES OUT ALL SOIL VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!

  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id
  integer :: NHW,NHE,NVN,NVS

  integer :: K,M,NX,NY,N,L,NGL,idsa
!     execution begins here
  integer, target  :: datic_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrc_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrc_2d(bounds%begc:bounds%endc,1:JZ+1)
  real(r8), target, allocatable :: datrc_3d(:,:,:)
  real(r8), target, allocatable :: datrc_4d(:,:,:,:)
  real(r8), target, allocatable :: datrc_5d(:,:,:,:,:)
  integer, pointer :: dat1pr(:)
  real(r8),pointer :: datpr1(:),datpr2(:,:),datpr3(:,:,:)
  real(r8),pointer :: datpr4(:,:,:,:),datpr5(:,:,:,:,:)
  integer :: sz3,sz4,sz5
  integer :: ncols, npfts  
  integer :: ndoms
  ncols = bounds%ncols
  npfts = bounds%npfts

  ndoms =trc_confs%NDOMS
  sz3=max(trc_confs%NGasTracers,  &
          trc_confs%NSolutTracers,  &
          trc_confs%NSaltTracers, & 
          trc_confs%NPrecipTracers,  &
          trc_confs%nxtracers,  & 
          trc_confs%NnutrientTracers,  &
          NumMicrobAutrophCmplx,jcplx)
  sz4=NumLiveHeterBioms
  sz5=jsken
  allocate(datrc_3d(bounds%begc:bounds%endc,1:sz3,1:JZ+1))
  allocate(datrc_4d(bounds%begc:bounds%endc,1:sz4,1:sz3,1:JZ+1))    
  allocate(datrc_5d(bounds%begc:bounds%endc,1:sz5,1:sz4,1:sz3,1:JZ+1))    

  NHW = bounds%NHW;NVN = bounds%NVN
  NHE = bounds%NHE;NVS = bounds%NVS

  call restartvar(ncid, flag, varname='CRAIN', &
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=CRAIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TSedmErossLoss_lnds', &
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=TSedmErossLoss_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TAREA', &
       long_name='total surface grid area', units='m2', &
       interpinic_flag='skip', data=TAREA, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HEATIN_lnd', &
       long_name='total surface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HEATIN_lnd,missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_O2_lnd', &
       long_name='total surface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_O2_lnd, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TORGF', &
       long_name='total organic C amendment', units='g d-2', &
       interpinic_flag='skip', data=TORGF, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TORGN', &
       long_name='total organic N amendment', units='g d-2', &
       interpinic_flag='skip', data=TORGN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TORGP', &
       long_name='total organic P amendment', units='g d-2', &
       interpinic_flag='skip', data=TORGP, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_N2_lnd', &
       long_name='total surface N2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_N2_lnd, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_CO2_lnd', &
       long_name='total surface CO2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_CO2_lnd, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='QH2OLoss_lnds', &
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=QH2OLoss_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CEVAP', &
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=CEVAP, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CRUN', &
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=CRUN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HeatOut_lnds', &
       long_name='total subsurface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HeatOut_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='OXYGOU', &
       long_name='total subsurface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=OXYGOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TCOU', &
       long_name='total subsurface C flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TZOU', &
       long_name='total subsurface N flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmn), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TPOU', &
       long_name='total subsurface P flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmp), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TZIN', &
       long_name='total surface N flux', units='g d-2', &
       interpinic_flag='skip', data=TZIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TPIN', &
       long_name='total surface P flux', units='g d-2', &
       interpinic_flag='skip', data=TPIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='XCSN', &
       long_name='total LitrFall C', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='XZSN', &
       long_name='total LitrFall N', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmn), missing_value=spval, &
       fill_value=spval)

  
  call restartvar(ncid, flag, varname='XPSN', &
       long_name='total LitrFall P', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmp), missing_value=spval, &
       fill_value=spval)
  
  
  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFNHB', dim1name='column',&
       long_name='banded NH4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IFNHB,datic_1d)
  else
    !print*,'IFNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IFNHB,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFNHB', dim1name='column',&
       long_name='banded NH4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval,flag_values=(/0,1/))

  endif

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NU', dim1name='column',&
       long_name='upper soil layer id', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NU,datic_1d)
  else 
    !print*,'NP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NU,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NU', dim1name='column',&
       long_name='upper soil layer id', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IDWaterTable', dim1name='column',&
       long_name='water table flag from site file', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IDWaterTable,datic_1d)
  else
    !print*,'IDWaterTable'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IDWaterTable,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IDWaterTable', dim1name='column',&
       long_name='water table flag from site file', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)

  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFNOB', dim1name='column',&
       long_name='banded NO3 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IFNOB,datic_1d)
  else
    !print*,'IFNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IFNOB,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFNOB', dim1name='column',&
       long_name='banded NO3 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFPOB', dim1name='column',&
       long_name='banded H2PO4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IFPOB,datic_1d)
  else
    !print*,'IFPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IFPOB,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IFPOB', dim1name='column',&
       long_name='banded H2PO4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IUTYP', dim1name='column',&
       long_name='urea hydrolysis inhibitor type (1=no,2=yes)', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IUTYP,datic_1d)
  else
    !print*,'IUTYP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IUTYP,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IUTYP', dim1name='column',&
       long_name='urea hydrolysis inhibitor type (1=no,2=yes)', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='CanopyHeight_col', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeight_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeight_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='CanopyHeight_col', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='ZERO4PlantDisplace_col', dim1name='column',&
       long_name='threshold for plant shoot displacement', units='', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,ZERO4PlantDisplace_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZERO4PlantDisplace_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='ZERO4PlantDisplace_col', dim1name='column',&
       long_name='threshold for plant shoot displacement', units='', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='AbvCanopyBndlResist_col', dim1name='column',&
       long_name='above canopy boundary layer resistance', units='h m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,AbvCanopyBndlResist_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AbvCanopyBndlResist_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='AbvCanopyBndlResist_col', dim1name='column',&
       long_name='above canopy boundary layer resistance', units='h m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='PrecIntceptByCanopy_col', dim1name='column',&
       long_name='net water transfer to whole grid canopy', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,PrecIntceptByCanopy_col,datrc_1d)
  else
    !print*,'PrecIntceptByCanopy_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PrecIntceptByCanopy_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='PrecIntceptByCanopy_col', dim1name='column',&
       long_name='net water transfer to whole grid canopy', units='MJ d-2 t-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='TSED', dim1name='column',&
       long_name='erosion rate', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,TSED,datrc_1d)
  else
    !print*,'TSED'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TSED,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='TSED', dim1name='column',&
       long_name='erosion rate', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='SoilSurfRoughnesst0_col', dim1name='column',&
       long_name='initial soil surface roughness height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilSurfRoughnesst0_col,datrc_1d)
  else
    !print*,'SoilSurfRoughnesst0_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilSurfRoughnesst0_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SoilSurfRoughnesst0_col', dim1name='column',&
       long_name='initial soil surface roughness height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='LWRadCanG', dim1name='column',&
       long_name='grid total canopy LW emission', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LWRadCanG,datrc_1d)
  else
    !print*,'LWRadCanG'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LWRadCanG,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='LWRadCanG', dim1name='column',&
       long_name='grid total canopy LW emission', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanWat_col', dim1name='column',&
       long_name='total water stored in dry matter', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanWat_col,datrc_1d)
  else
    !print*,'LWRadCanG'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanWat_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanWat_col', dim1name='column',&
       long_name='total water stored in dry matter', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyHeatStor_col', dim1name='column',&
       long_name='total canopy heat content', units='MJ d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeatStor_col,datrc_1d)
  else
    !print*,'LWRadCanG'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeatStor_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyHeatStor_col', dim1name='column',&
       long_name='total canopy heat content', units='MJ d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='LWRadBySurf_col', dim1name='column',&
       long_name='longwave radiation emitted from ground surface (including snow and litter)', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf_col,datrc_1d)
  else
    !print*,'LWRadBySurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='LWRadBySurf_col', dim1name='column',&
       long_name='longwave radiation emitted from ground surface  (including snow and litter)', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='Eco_NEE_col', dim1name='column',&
       long_name='total canopy net CO2 exchange', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NEE_col,datrc_1d)
  else
    !print*,'Eco_NEE_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NEE_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='Eco_NEE_col', dim1name='column',&
       long_name='total canopy net CO2 exchange', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='RIB', dim1name='column',&
       long_name='Richardson number', units=' ', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RIB,datrc_1d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIB,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='RIB', dim1name='column',&
       long_name='Richardson number', units=' ', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  


  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanH2OHeldVg_col', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CanH2OHeldVg_col,datrc_1d) 
  else
    !print*,'CanH2OHeldVg_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanH2OHeldVg_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanH2OHeldVg_col', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='VWatStoreCapSurf_col', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf_col,datrc_1d) 
  else
    !print*,'VWatStoreCapSurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf_col,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='VWatStoreCapSurf_col', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='QRain_CumYr_col', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QRain_CumYr_col,datrc_1d) 
  else
    !print*,'QRain_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QRain_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='QRain_CumYr_col', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyLeafArea_col', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_col,datrc_1d) 
  else
    !print*,'CanopyLeafArea_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyLeafArea_col', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='StemArea_col', dim1name='column',&
       long_name='grid canopy stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_col,datrc_1d) 
  else
    !print*,'StemArea_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='StemArea_col', dim1name='column',&
       long_name='total canopy stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='PPT', dim1name='column',&
       long_name='total plant population', units='# d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PlantPopu_col,datrc_1d) 
  else
    !print*,'PPT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PlantPopu_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='PPT', dim1name='column',&
       long_name='total plant population', units='# d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VLWatheldCapSurf_col', dim1name='column',&
       long_name='soil surface water retention capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatheldCapSurf_col,datrc_1d) 
  else
    !print*,'VLWatheldCapSurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatheldCapSurf_col,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='VLWatheldCapSurf_col', dim1name='column',&
       long_name='soil surface water retention capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SoiSurfRoughness', dim1name='column',&
       long_name='soil surface roughness height for calculating runoff velocity', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughness,datrc_1d) 
  else
    !print*,'SoiSurfRoughness'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughness,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='SoiSurfRoughness', dim1name='column',&
       long_name='soil surface roughness height for calculating runoff velocity', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)    
    call restartvar(ncid, flag, varname='SurfGasFlx_col', dim1name='column',&
       dim2name='gastrcs',long_name='total soil gas flux', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SurfGasFlx_col,datrc_2d) 
  else
    !print*,'UCH4G'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SurfGasFlx_col,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)          
    call restartvar(ncid, flag, varname='SurfGasFlx_col', dim1name='column',&
       dim2name='gastrcs',long_name='total soil gas flux', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CO2byFire_CumYr_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_CumYr_col,datrc_1d) 
  else
    !print*,'CO2byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_CumYr_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CO2byFire_CumYr_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    

  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CH4byFire_CumYr_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_CumYr_col,datrc_1d) 
  else
    !print*,'CH4byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CH4byFire_CumYr_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='O2byFire_CumYr_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_CumYr_col,datrc_1d) 
  else
    !print*,'O2byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='O2byFire_CumYr_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_CumYr_col,datrc_1d) 
  else
    !print*,'N2ObyFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_CumYr_col,datrc_1d)
    datpr1 => datrc_1d             
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_CumYr_col,datrc_1d) 
  else
    !print*,'NH3byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_CumYr_col,datrc_1d) 
  else
    !print*,'PO4byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='AmendCFlx_CumYr_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,AmendCFlx_CumYr_col,datrc_1d) 
  else
    !print*,'AmendCFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AmendCFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='AmendCFlx_CumYr_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FertNFlx_CumYr_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FertNFlx_CumYr_col,datrc_1d) 
  else
    !print*,'FertNFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertNFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FertNFlx_CumYr_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FerPFlx_CumYr_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FerPFlx_CumYr_col,datrc_1d) 
  else
    !print*,'FerPFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FerPFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FerPFlx_CumYr_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='H2OLoss_CumYr_col', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,H2OLoss_CumYr_col,datrc_1d) 
  else
    !print*,'H2OLoss_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,H2OLoss_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='H2OLoss_CumYr_col', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='QEvap_CumYr_col', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QEvap_CumYr_col,datrc_1d) 
  else
    !print*,'QEvap_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QEvap_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='QEvap_CumYr_col', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='Qrunoff_CumYr_col', dim1name='column',&
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Qrunoff_CumYr_col,datrc_1d) 
  else
    !print*,'Qrunoff_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Qrunoff_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='Qrunoff_CumYr_col', dim1name='column',&
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDOCFlx_col', dim1name='column',&
       long_name='total surface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOCFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDOCFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOCFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDOCFlx_col', dim1name='column',&
       long_name='total surface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDOCFlx_col', dim1name='column',&
       long_name='total subsurface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOCFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDOCFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOCFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDOCFlx_col', dim1name='column',&
       long_name='total subsurface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDONFlx_CumYr_col', dim1name='column',&
       long_name='total surface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDONFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_CumYr_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDONFlx_CumYr_col', dim1name='column',&
       long_name='total surface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    

  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDONFlx_col', dim1name='column',&
       long_name='total subsurface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDONFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDONFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDONFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDONFlx_col', dim1name='column',&
       long_name='total subsurface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDOPFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDOPFlx_col', dim1name='column',&
       long_name='total subsurface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOPFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDOPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOPFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDOPFlx_col', dim1name='column',&
       long_name='total subsurface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDICFlx_col', dim1name='column',&
       long_name='total surface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDICFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDICFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDICFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDICFlx_col', dim1name='column',&
       long_name='total surface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDICFlx_col', dim1name='column',&
       long_name='total subsurface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDICFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDICFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDICFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDICFlx_col', dim1name='column',&
       long_name='total subsurface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSufDINFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDINFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDINFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDINFlx_col', dim1name='column',&
       long_name='total subsurface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDINFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDINFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDINFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDINFlx_col', dim1name='column',&
       long_name='total subsurface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDIPFlx_CumYr_col'
   if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDIPFlx_col', dim1name='column',&
       long_name='total subsurface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDIPFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDIPFlx_col'
   if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDIPFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDIPFlx_col', dim1name='column',&
       long_name='total subsurface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)            
    call restartvar(ncid, flag, varname='LiterfalOrgM_col', dim1name='column',&
       dim2name='elmnts',long_name='total litrFall organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col(:,1:JY0,1:JX0),datrc_2d) 
  else
!    print*,'LiterfalOrgC_col',size(LiterfalOrgM_col,1),size(LiterfalOrgM_col,2)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col(:,1:JY0,1:JX0),datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LiterfalOrgM_col', dim1name='column',&
       dim2name='elmnts',long_name='total litrFall organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)            
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_col', dim1name='column',&
       dim2name='elmnts',long_name='total harvested organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,EcoHavstElmnt_CumYr_col(:,1:JY0,1:JX0),datrc_2d) 
  else
!    print*,'LiterfalOrgC_col',size(LiterfalOrgM_col,1),size(LiterfalOrgM_col,2)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EcoHavstElmnt_CumYr_col(:,1:JY0,1:JX0),datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_col', dim1name='column',&
       dim2name='elmnts',long_name='total harvested organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='RootResp_CumYr_col', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RootResp_CumYr_col,datrc_1d) 
  else
    !print*,'RootResp_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RootResp_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='RootResp_CumYr_col', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='QDrain_col', dim1name='column',&
       long_name='total water drainage below root zone', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QDrain_col,datrc_1d) 
  else
    !print*,'QDrain_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QDrain_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='QDrain_col', dim1name='column',&
       long_name='total water drainage below root zone', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ZDRAIN', dim1name='column',&
       long_name='total N drainage below root zone', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ZDRAIN,datrc_1d) 
  else
    !print*,'ZDRAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZDRAIN,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ZDRAIN', dim1name='column',&
       long_name='total N drainage below root zone', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='PDRAIN', dim1name='column',&
       long_name='total P drainage below root zone', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PDRAIN,datrc_1d) 
  else
    !print*,'PDRAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PDRAIN,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='PDRAIN', dim1name='column',&
       long_name='total P drainage below root zone', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPNH4', dim1name='column',&
       long_name='total depth of NH4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNH4_col,datrc_1d) 
  else
    !print*,'DPNH4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNH4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPNH4', dim1name='column',&
       long_name='total depth of NH4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPNO3', dim1name='column',&
       long_name='total depth of NO3 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNO3_col,datrc_1d) 
  else
    !print*,'DPNO3'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNO3_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPNO3', dim1name='column',&
       long_name='total depth of NO3 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPPO4', dim1name='column',&
       long_name='total depth of PO4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthPO4_col,datrc_1d) 
  else
    !print*,'DPPO4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthPO4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPPO4', dim1name='column',&
       long_name='total depth of PO4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion1', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion1,datrc_1d) 
  else
    !print*,'SoilDetachability4Erosion1'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion1,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion1', dim1name='column',&
       long_name='soil detachability', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion2', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion2,datrc_1d) 
  else
    !print*,'SoilDetachability4Erosion2'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion2,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion2', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='CER', dim1name='column',&
       long_name='microbial C  erosion', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CER,datrc_1d) 
  else
    !print*,'CER'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CER,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='CER', dim1name='column',&
       long_name='microbial C  erosion', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='XER', dim1name='column',&
       long_name='soil detachment/deposition', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,XER,datrc_1d) 
  else
    !print*,'XER'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,XER,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='XER', dim1name='column',&
       long_name='soil detachment/deposition', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SedmErossLoss_CumYr_col', dim1name='column',&
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SedmErossLoss_CumYr_col,datrc_1d) 
  else
    !print*,'SedmErossLoss_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SedmErossLoss_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SedmErossLoss_CumYr_col', dim1name='column',&
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWN', dim1name='column',&
       long_name='row spacing of NH4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWN,datrc_1d) 
  else
    !print*,'ROWN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWN,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWN', dim1name='column',&
       long_name='row spacing of NH4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWO', dim1name='column',&
       long_name='row spacing of NO3 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWO,datrc_1d) 
  else
    !print*,'ROWO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWO,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWO', dim1name='column',&
       long_name='row spacing of NO3 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWP', dim1name='column',&
       long_name='row spacing of PO4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWP,datrc_1d) 
  else
    !print*,'ROWP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWP,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWP', dim1name='column',&
       long_name='row spacing of PO4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ExtWaterTablet0', dim1name='column',&
       long_name='initial external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTablet0,datrc_1d) 
  else
    !print*,'ExtWaterTablet0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTablet0,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ExtWaterTablet0', dim1name='column',&
       long_name='initial external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DTBLD', dim1name='column',&
       long_name='depth of artificial water table adjusted for elevation', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DTBLD,datrc_1d) 
  else
    !print*,'DTBLD'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DTBLD,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DTBLD', dim1name='column',&
       long_name='depth of artificial water table adjusted for elevation', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_NBP_CumYr_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_NBP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NBP_CumYr_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='VLitR_col', dim1name='column',&
       long_name='surface litter volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLitR_col,datrc_1d) 
  else
    !print*,'VLitR_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLitR_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='VLitR_col', dim1name='column',&
       long_name='surface litter volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SED', dim1name='column',&
       long_name='sediment transport', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SED,datrc_1d) 
  else
    !print*,'SED'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SED,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SED', dim1name='column',&
       long_name='sediment transport', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_GPP_CumYr_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_GPP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_GPP_CumYr_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_AutoR_CumYr_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_AutoR_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_AutoR_CumYr_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_NPP_CumYr_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_NPP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NPP_CumYr_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_HR_CumYr_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_HR_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_HR_CumYr_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Air_Heat_Latent_store_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Latent_store_col,datrc_1d) 
  else
    !print*,'Air_Heat_Latent_store_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Latent_store_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Air_Heat_Latent_store_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Air_Heat_Sens_store_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Sens_store_col,datrc_1d) 
  else
    !print*,'Air_Heat_Sens_store_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Sens_store_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Air_Heat_Sens_store_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DayLensCurr_col', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLensCurr_col,datrc_1d) 
  else
    !print*,'DayLensCurr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLensCurr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DayLensCurr_col', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='DayLenthPrev_col', dim1name='column',&
       long_name='daylength of previous day', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev_col,datrc_1d) 
  else
    !print*,'DayLenthPrev_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev_col,datrc_1d)   
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='DayLenthPrev_col', dim1name='column',&
       long_name='daylength of previous day', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='ExtWaterTable', dim1name='column',&
       long_name='current external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable_col,datrc_1d) 
  else
    !print*,'ExtWaterTable'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable_col,datrc_1d)   
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='ExtWaterTable', dim1name='column',&
       long_name='current external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='DTBLY', dim1name='column',&
       long_name='artificial water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DTBLY,datrc_1d) 
  else
    !print*,'DTBLY'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DTBLY,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='DTBLY', dim1name='column',&
       long_name='artificial water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:micpar%NumOfLitrCmplxs)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nlitromcomplx',&
       long_name='surface OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RC0(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'RC0',size(RC0,1),size(RC0,2),jcplx
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RC0(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:micpar%NumOfLitrCmplxs)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nlitromcomplx',&
       long_name='surface OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE_col', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE_col,datrc_1d) 
  else
    !print*,'VcumDrySnoWE_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE_col', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow_col', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow_col,datrc_1d) 
  else
    !print*,'VcumWatSnow_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow_col', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow_col', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow_col,datrc_1d) 
  else
    !print*,'VcumIceSnow_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow_col', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI_col', dim1name='column',&
       long_name='total snowpack volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI_col,datrc_1d) 
  else
    !print*,'VcumSnoDWI_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI_col', dim1name='column',&
       long_name='total snowpack volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SnowDepth', dim1name='column',&
       long_name='snowpack depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth_col,datrc_1d) 
  else
    !print*,'SnowDepth'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SnowDepth', dim1name='column',&
       long_name='snowpack depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='EnergyImpact4Erosion', dim1name='column',&
       long_name='cumulative rainfall energy impact on soil surface', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,EnergyImpact4Erosion,datrc_1d) 
  else
    !print*,'EnergyImpact4Erosion'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EnergyImpact4Erosion,datrc_1d)    
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='EnergyImpact4Erosion', dim1name='column',&
       long_name='cumulative rainfall energy impact on soil surface', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
   endif 

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLDrySnoWE', dim1name='column',dim2name='levsno',&
       long_name='water equivalent dry snow in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLDrySnoWE'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)    
    call restartvar(ncid, flag, varname='VLDrySnoWE', dim1name='column',dim2name='levsno',&
       long_name='water equivalent dry snow in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLIceSnow', dim1name='column',dim2name='levsno',&
       long_name='snow ice volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLIceSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)    
    call restartvar(ncid, flag, varname='VLIceSnow', dim1name='column',dim2name='levsno',&
       long_name='snow ice volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLWatSnow', dim1name='column',dim2name='levsno',&
       long_name='snow water volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLWatSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLWatSnow', dim1name='column',dim2name='levsno',&
       long_name='snow water volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWIprev_snvr', dim1name='column',dim2name='levsno',&
       long_name='total snow volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWIprev_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLSnoDWIprev_snvr
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWIprev_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWIprev_snvr', dim1name='column',dim2name='levsno',&
       long_name='total snow volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnoDensL', dim1name='column',dim2name='levsno',&
       long_name='snowpack density', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnoDens_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'SnoDensL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnoDens_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnoDensL', dim1name='column',dim2name='levsno',&
       long_name='snowpack density', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowThickL_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowThickL_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'SnowThickL_snvr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowThickL_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowThickL_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLHeatCapSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TKSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'TKSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TKSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TCSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'TCSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TCSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JS)    
    call restartvar(ncid, flag, varname='trcg_solsml', dim1name='column',dim2name='gastrcs',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml_snvr(:,:,1:JY0,1:JX0),datrc_3d) 
  else
    !print*,'trcg_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml_snvr(:,:,1:JY0,1:JX0),datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JS)        
    call restartvar(ncid, flag, varname='trcg_solsml', dim1name='column',dim2name='gastrcs',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif 

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NnutrientTracers,1:JS)            
    call restartvar(ncid, flag, varname='trcn_solsml', dim1name='column',dim2name='NnutrientTracers',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml_snvr,datrc_3d) 
  else
    !print*,'trcn_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml_snvr,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NnutrientTracers,1:JS)            
    call restartvar(ncid, flag, varname='trcn_solsml', dim1name='column',dim2name='NnutrientTracers',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SoilFracAsMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='macropore fraction', units='m3/m3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP_vr,datrc_2d) 
  else
    !print*,'SoilFracAsMacP_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoilFracAsMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='thickness of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:3,1:JZ+1)    
    call restartvar(ncid, flag, varname='DLYR', dim1name='column',dim2name='sdim',&
       dim3name='levsoi1',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DLYR,datrc_3d) 
  else
    !print*,'DLYR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DLYR,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:3,1:JZ+1)            
    call restartvar(ncid, flag, varname='DLYR', dim1name='column',dim2name='sdim',&
       dim3name='levsoi1',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='CumDepz2LayerBot_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumDepz2LayerBot_vr,datrc_2d) 
  else
    !print*,'CumDepz2LayerBot_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumDepz2LayerBot_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumDepz2LayerBot_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness_vr,datrc_2d) 
  else
    !print*,'CumSoilThickness_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoiBulkDensityt0_vr', dim1name='column',dim2name='levsoi',&
       long_name='initial bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0_vr,datrc_2d) 
  else
    !print*,'SoiBulkDensityt0_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='SoiBulkDensityt0_vr', dim1name='column',dim2name='levsoi',&
       long_name='initial bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoiBulkDensity', dim1name='column',dim2name='levsoi1',&
       long_name='soil bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensity_vr,datrc_2d) 
  else
    !print*,'SoiBulkDensity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensity_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoiBulkDensity', dim1name='column',dim2name='levsoi1',&
       long_name='soil bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='CSoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
       long_name='soil organic matter content density', units='g kg-1 m-3', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CSoilOrgM_vr,datrc_3d) 
  else
    !print*,'CORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CSoilOrgM_vr,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='CSoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
       long_name='soil organic matter content density', units='g kg-1 m-3', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,POROS_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,POROS_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PH', dim1name='column',dim2name='levsoi1',&
       long_name='soil pH', units='none', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PH,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PH,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PH', dim1name='column',dim2name='levsoi1',&
       long_name='soil pH', units='none', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicPMass_vr', dim1name='column',dim2name='levsoi1',&
       long_name='mass of soil layer', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicPMass_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicPMass_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicPMass_vr', dim1name='column',dim2name='levsoi1',&
       long_name='mass of soil layer', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilPoreMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilPoreMicP_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilPoreMicP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilPoreMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PSISoilMatricP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='MPa', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PSISoilMatricP_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PSISoilMatricP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PSISoilMatricP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil matric potential', units='MPa', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETW_vr', dim1name='column',dim2name='levsoi1',&
       long_name='volumetric water content', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,THETW_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,THETW_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETW_vr', dim1name='column',dim2name='levsoi1',&
       long_name='volumetric water content', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity_vr,datrc_2d) 
  else
    !print*,'FieldCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint_vr', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint_vr,datrc_2d)   
  else
    !print*,'WiltPoint_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint_vr,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint_vr', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert_vr,datrc_2d)   
  else
    !print*,'SatHydroCondVert_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil horizontal saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn_vr,datrc_2d)   
  else
    !print*,'SatHydroCondHrzn_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil horizontal saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SAND', dim1name='column',dim2name='levsoi',&
       long_name='soil sand content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SAND,datrc_2d)   
  else
    !print*,'SAND'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SAND,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SAND', dim1name='column',dim2name='levsoi',&
       long_name='soil sand content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SILT', dim1name='column',dim2name='levsoi',&
       long_name='soil silt content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SILT,datrc_2d)   
  else
    !print*,'SILT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SILT,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SILT', dim1name='column',dim2name='levsoi',&
       long_name='soil silt content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CLAY', dim1name='column',dim2name='levsoi',&
       long_name='soil clay content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,CLAY,datrc_2d)   
  else
    !print*,'CLAY'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CLAY,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CLAY', dim1name='column',dim2name='levsoi',&
       long_name='soil clay content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP_vr,datrc_2d)   
  else
    !print*,'VLWatMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicPX', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content before wetting front', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX_vr,datrc_2d)     
  else
    !print*,'VLWatMicPX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicPX', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content before wetting front', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLiceMicP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP_vr,datrc_2d)     
  else
    !print*,'VLiceMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLiceMicP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLsoiAirP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore air content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP_vr,datrc_2d)     
  else
    !print*,'VLsoiAirP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLsoiAirP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore air content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP_vr,datrc_2d)     
  else
    !print*,'VLMicP_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP_vr,datrc_2d)     
  else
    !print*,'VLSoilMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETY_vr', dim1name='column',dim2name='levsoi1',&
       long_name='air-dry water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,THETY_vr,datrc_2d)     
  else
    !print*,'VLSoilMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,THETY_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETY_vr', dim1name='column',dim2name='levsoi1',&
       long_name='air-dry water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)    
    call restartvar(ncid, flag, varname='trcx_solml_vr', dim1name='column',dim2name='xtracers',&
       dim3name='levsoi1',long_name='exchangeable tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml_vr,datrc_3d)     
  else
    !print*,'trcx_solml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)        
    call restartvar(ncid, flag, varname='trcx_solml_vr', dim1name='column',dim2name='xtracers',&
       dim3name='levsoi1',long_name='exchangeable tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='TCS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,TCS,datrc_2d)     
  else
    !print*,'TCS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TCS,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TCS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TKS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,TKS_vr,datrc_2d)     
  else
    !print*,'TKS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKS_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TKS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  


  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ d-2 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity_vr,datrc_2d)     
  else
    !print*,'VHeatCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ d-2 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacitySoilM_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil solid heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacitySoilM_vr,datrc_2d)     
  else
    !print*,'VHeatCapacitySoilM_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacitySoilM_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacitySoilM_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil solid heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)    
    call restartvar(ncid, flag, varname='trc_gasml_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trc_gasml_vr,datrc_3d)     
  else
    !print*,'trc_gasml_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trc_gasml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)        
    call restartvar(ncid, flag, varname='trc_gasml_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trc_solml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trc_solml_vr,datrc_3d)     
  else
    !print*,'trc_solml_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trc_solml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trc_solml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)           
    call restartvar(ncid, flag, varname='trc_soHml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trc_soHml_vr,datrc_3d)     
  else
    !print*,'trc_soHml_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trc_soHml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trc_soHml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   

  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2GasXchangePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2GasXchangePrev_vr,datrc_2d)     
  else
    !print*,'RO2GasXchangePrev_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2GasXchangePrev_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2GasXchangePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCO2GasFlxPrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CO2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RCO2GasFlxPrev_vr,datrc_2d)     
  else
    !print*,'RCO2GasFlxPrev_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCO2GasFlxPrev_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCO2GasFlxPrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CO2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2AquaSourcePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2AquaSourcePrev_vr,datrc_2d)       
  else
    !print*,'RO2AquaSourcePrev_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2AquaSourcePrev_vr,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2AquaSourcePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4F', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RCH4F,datrc_2d)       
  else
    !print*,'RCH4F'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4F,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4F', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4L', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RCH4PhysexchPrev_vr,datrc_2d)       
  else
    !print*,'RCH4L'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4PhysexchPrev_vr,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4L', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoO2DmndResp_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoO2DmndResp_vr,datrc_2d)      
  else
    !print*,'REcoO2DmndResp_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoO2DmndResp_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoO2DmndResp_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoNH4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoNO3DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkSoil_vr,datrc_2d)      
  else
    !print*,'RNO2EcoUptkSoil_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OEcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2OEcoUptkSoil_vr,datrc_2d)      
  else
    !print*,'RN2OEcoUptkSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2OEcoUptkSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OEcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoH2PO4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoH1PO4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoNH4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoNO3DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkBand_vr,datrc_2d)      
  else
    !print*,'RNO2EcoUptkBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoH2PO4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoH1PO4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)        
    call restartvar(ncid, flag, varname='RDOMEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RDOMEcoDmndK_vr,datrc_3d)      
  else
    !print*,'RDOMEcoDmndK_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RDOMEcoDmndK_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RDOMEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)        
    call restartvar(ncid, flag, varname='FracBulkSOMC_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='mass fraction of complex OM in bulk OM', units='none', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,FracBulkSOMC_vr,datrc_3d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FracBulkSOMC_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='FracBulkSOMC_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='mass fraction of complex OM in bulk OM', units='none', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RAcetateEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial acetate uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateEcoDmndK_vr,datrc_3d)      
  else
    !print*,'RAcetateEcoDmndK_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateEcoDmndK_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RAcetateEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial acetate uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VLWatMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP_vr,datrc_2d)      
  else
    !print*,'VLWatMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLWatMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         

  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLiceMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP_vr,datrc_2d)      
  else
    !print*,'VLiceMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLiceMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                
    call restartvar(ncid, flag, varname='VLMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMacP_vr,datrc_2d)      
  else
    !print*,'VLMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMacP_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='VLMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='totRootLenDens_vr', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,totRootLenDens_vr,datrc_2d)      
  else
    !print*,'totRootLenDens_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,totRootLenDens_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='totRootLenDens_vr', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                    
    call restartvar(ncid, flag, varname='CanopyHeightZ_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightZ_col,datrc_2d)      
  else
    !print*,'ZL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='CanopyHeightZ_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)                        
    call restartvar(ncid, flag, varname='CanopyLeafAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafAareZ_col,datrc_2d)      
  else
    !print*,'CanopyLeafAareZ_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafAareZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyLeafAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='TAU_DirRadTransm', dim1name='column',dim2name='levcan1',&
       long_name='Direct radiation transmission coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirRadTransm,datrc_2d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirRadTransm,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)    
    call restartvar(ncid, flag, varname='TAU_DirRadTransm', dim1name='column',dim2name='levcan1',&
       long_name='Direct radiation transmission coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='TAU_RadThru', dim1name='column',dim2name='levcan1',&
       long_name='Through canopy radiation coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,TAU_RadThru,datrc_2d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TAU_RadThru,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)    
    call restartvar(ncid, flag, varname='TAU_RadThru', dim1name='column',dim2name='levcan1',&
       long_name='Through canopy radiation coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemAareZ_col,datrc_2d)      
  else
    !print*,'CanopyStemAareZ_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemAareZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='tCanLeafC_cl', dim1name='column',dim2name='levcan',&
       long_name='total leaf mass', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,tCanLeafC_cl,datrc_2d)      
  else
    !print*,'tCanLeafC_cl'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,tCanLeafC_cl,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='tCanLeafC_cl', dim1name='column',dim2name='levcan',&
       long_name='total leaf mass', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcs_VLN_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='nutrient salt tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,trcs_VLN_vr,datrc_3d)      
  else
    !print*,'trcs_VLN_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_VLN_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcs_VLN_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='nutrient salt tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NPrecipTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcp_saltpml_vr', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,trcp_saltpml_vr,datrc_3d)      
  else
    !print*,'trcp_salml'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trcp_saltpml_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NPrecipTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trcp_saltpml_vr', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)        
    call restartvar(ncid, flag, varname='SolutesIonConc_vr', dim1name='column',dim2name='levcan',&
       long_name='solution ion concentratiom', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,SolutesIonConc_vr,datrc_2d)      
  else
    !print*,'SolutesIonConc_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolutesIonConc_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)           
    call restartvar(ncid, flag, varname='SolutesIonConc_vr', dim1name='column',dim2name='levcan',&
       long_name='solution ion concentratiom', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  


  IF(salt_model)THEN

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)           
      call restartvar(ncid, flag, varname='trcSalt_solml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi1',long_name='soil aqueous salt content micropre', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml_vr,datrc_3d)      
    else
      !print*,'trcSalt_solml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml_vr,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)               
      call restartvar(ncid, flag, varname='trcSalt_solml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi1',long_name='soil aqueous salt content micropre', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JS)               
      call restartvar(ncid, flag, varname='trcs_solsml', dim1name='column',dim2name='satracers',&
         dim3name='levsno',long_name='snowpack salt dissolved tracers', units='mol d-2', &
         interpinic_flag='skip', data=datpr3, missing_value=spval, &
         fill_value=spval)  
        call cpcol(flag,NHW,NHE,NVN,NVS,trc_Saltml_snvr,datrc_3d)      
    else
      !print*,'trcs_solsml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trc_Saltml_snvr,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JS)    
      call restartvar(ncid, flag, varname='trcs_solsml', dim1name='column',dim2name='satracers',&
         dim3name='levsno',long_name='snowpack salt dissolved tracers', units='mol d-2', &
         interpinic_flag='skip', data=datpr3, missing_value=spval, &
         fill_value=spval)  

    endif  

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml_vr,datrc_3d)      
    else
      !print*,'trcSalt_soHml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml_vr,datrc_3d)          
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
    endif  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2DmndHetert', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndHetert,datrc_4d)      
  else
    !print*,'RO2DmndHetert'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndHetert,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2DmndHetert', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO3ReduxDmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO2DmndReduxSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total heterotrophic microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxHeter_vr,datrc_4d)      
  else
    !print*,'RN2ODmndReduxHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total heterotrophic microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO3ReduxDmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO2DmndReduxBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNH4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO3DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNH4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO3DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H1PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H1PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RDOCUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RDOCUptkHeter_vr,datrc_4d)      
  else
    !print*,'RDOCUptkHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RDOCUptkHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RDOCUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RAcetateUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateUptkHeter_vr,datrc_4d)      
  else
    !print*,'RAcetateUptkHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateUptkHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetrMicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RAcetateUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNH4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RNH4DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNH4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNO3DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RNO3DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNO3DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RH2PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RH2PO4DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetrMicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RH2PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='mBiomeHeter_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element',units='g d-2', interpinic_flag='skip', &
      data=datpr5, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeHeter_vr,datrc_5d)      
  else
    !print*,'OMC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeHeter_vr,datrc_5d)        
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='mBiomeHeter_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  
  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2DmndAutort_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndAutort_vr,datrc_3d)      
  else
    !print*,'RO2DmndAutort'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndAutort_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2DmndAutort_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutor', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutor,datrc_3d)      
  else
    !print*,'RNH3OxidAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutor,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutor', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2OxidAutor', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2OxidAutor,datrc_3d)      
  else
    !print*,'RNO2OxidAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2OxidAutor,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2OxidAutor', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxAutor_vr,datrc_3d)      
  else
    !print*,'RN2ODmndReduxAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutorBand', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutorBand,datrc_3d)      
  else
    !print*,'RNH3OxidAutorBand'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutorBand,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutorBand', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2OxidAutorBand', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2OxidAutorBand,datrc_3d)     
  else
    !print*,'RNO2OxidAutorBand'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2OxidAutorBand,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2OxidAutorBand', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RNH4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RNO3UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RNH4UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RNO3UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)    
    call restartvar(ncid, flag, varname='RNH4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RNH4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RNH4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RNO3UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RNO3UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RNO3UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RH2PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H2PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RH2PO4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RH2PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H2PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RH1PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H1PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RH1PO4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutrophCmplx)       
    call restartvar(ncid, flag, varname='RH1PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H1PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)       
    call restartvar(ncid, flag, varname='mBiomeAutor_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeAutor_vr,datrc_4d)     
  else
    !print*,'OMCff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeAutor_vr,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)           
    call restartvar(ncid, flag, varname='mBiomeAutor_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='OMBioResdu_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OMBioResdu_vr,datrc_5d)     
  else
    !print*,'ORC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMBioResdu_vr,datrc_5d)       
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='OMBioResdu_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)               
    call restartvar(ncid, flag, varname='DOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic C micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_vr,datrc_4d)     
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_vr,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      

  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MacP', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter macropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP_vr,datrc_4d)
  else
    !print*,'OQCH'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MacP', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter macropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SorbedOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SorbedOM_vr,datrc_4d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SorbedOM_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SorbedOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      

  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SolidOM_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SolidOM_vr,datrc_5d)
  else
    !print*,'OSC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolidOM_vr,datrc_5d)  
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOM_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOMAct_vr', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SolidOMAct_vr,datrc_4d)
  else
    !print*,'OSA'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolidOMAct_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOMAct_vr', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndSoilChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndSoilChemo_vr,datrc_2d)
  else
    !print*,'RNO2DmndSoilChemo_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndSoilChemo_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndSoilChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndBandChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndBandChemo_vr,datrc_2d)
  else
    !print*,'RNO2DmndSoilChemo_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndBandChemo_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndBandChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  
  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
      long_name='total soil organic matter', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      

    call cpcol(flag,NHW,NHE,NVN,NVS,SoilOrgM_vr,datrc_3d)

  else
    !print*,'ORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilOrgM_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)   
    call restartvar(ncid, flag, varname='SoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
      long_name='total soil organic matter', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMLitrC_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OMLitrC_vr,datrc_2d)
  else
    !print*,'ORGR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMLitrC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMLitrC_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_soil_vr', dim1name='column',dim2name='fertN',&
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_soil_vr,datrc_3d)
  else
    !print*,'FertN_soil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_soil_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_soil_vr', dim1name='column',dim2name='fertN', &
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_Band_vr', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_Band_vr,datrc_3d)
  else
    !print*,'FertN_Band_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_Band_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_Band_vr', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNH4_vr,datrc_2d)
  else
    !print*,'WDNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNH4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='AEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='anion exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,AEC_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AEC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='AEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='anion exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='cation exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CEC_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CEC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='cation exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CFE', dim1name='column', &
      dim2name='levsoi',long_name='soil Fe content', units='mg Fe kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CFE,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CFE,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CFE', dim1name='column', &
      dim2name='levsoi',long_name='soil Fe content', units='mg Fe kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCA', dim1name='column', &
      dim2name='levsoi',long_name='soil Ca content', units='mg Ca kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CCA,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CCA,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCA', dim1name='column', &
      dim2name='levsoi',long_name='soil Ca content', units='mg Ca kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CAL', dim1name='column', &
      dim2name='levsoi',long_name='soil Al content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CAL,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CAL,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CAL', dim1name='column', &
      dim2name='levsoi',long_name='soil Al content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNHB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNH4_vr,datrc_2d)
  else
    !print*,'DPNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNH4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNHB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNOB', dim1name='column', &
      dim2name='levsoi',long_name='width of NO3 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNO3_vr,datrc_2d)
  else
    !print*,'WDNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNO3_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNOB', dim1name='column', &
      dim2name='levsoi',long_name='width of NO3 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNO3_vr,datrc_2d)
  else
    !print*,'DPNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNO3_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDPOB', dim1name='column', &
      dim2name='levsoi',long_name='width of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthPO4_vr,datrc_2d)
  else
    !print*,'WDPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthPO4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDPOB', dim1name='column', &
      dim2name='levsoi',long_name='width of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPPOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessPO4_vr,datrc_2d)
  else
     !print*,'DPPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessPO4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPPOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHUI', dim1name='column', &
      dim2name='levsoi1',long_name='current inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNHUI,datrc_2d)
  else
    !print*,'ZNHUI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNHUI,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHUI', dim1name='column', &
      dim2name='levsoi1',long_name='current inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHU0', dim1name='column', &
      dim2name='levsoi1',long_name='urea hydrolysis inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNHU0,datrc_2d)
  else
    !print*,'ZNHU0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNHU0,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHU0', dim1name='column', &
      dim2name='levsoi1',long_name='urea hydrolysis inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFNI', dim1name='column', &
      dim2name='levsoi1',long_name='current nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNFNI,datrc_2d)  
  else
    !print*,'ZNFNI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNFNI,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFNI', dim1name='column', &
      dim2name='levsoi1',long_name='current nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFN0', dim1name='column', &
      dim2name='levsoi1',long_name='initial nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNFN0,datrc_2d)  
  else
    !print*,'ZNFN0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNFN0,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFN0', dim1name='column', &
      dim2name='levsoi1',long_name='initial nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPX', dim1name='column', &
      dim2name='month',long_name='accumulated change  for maximum temperature', &
      units='K', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDTPX,datrc_2d)    
  else
    !print*,'TDTPX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDTPX,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPX', dim1name='column', &
      dim2name='month',long_name='accumulated change for maximum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDRAD', dim1name='column', &
      dim2name='month',long_name='accumulated change for radiation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDRAD,datrc_2d)    
  else
    !print*,'TDRAD'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDRAD,datrc_2d)      
    datpr2 => datrc_2d (1:ncols,1:12)   
    call restartvar(ncid, flag, varname='TDRAD', dim1name='column', &
      dim2name='month',long_name='accumulated change for radiation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPN', dim1name='column', &
      dim2name='month',long_name='accumulated change for minimum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDTPN,datrc_2d)    
  else
    !print*,'TDTPN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDTPN,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPN', dim1name='column', &
      dim2name='month',long_name='accumulated change for minimum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDWND', dim1name='column', &
      dim2name='month',long_name='accumulated change for wind speed', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDWND,datrc_2d)    
  else
    !print*,'TDWND'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDWND,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)   
    call restartvar(ncid, flag, varname='TDWND', dim1name='column', &
      dim2name='month',long_name='accumulated change for wind speed', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDHUM', dim1name='column', &
      dim2name='month',long_name='accumulated change for humidity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDHUM,datrc_2d)    
  else
    !print*,'TDHUM'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDHUM,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDHUM', dim1name='column', &
      dim2name='month',long_name='accumulated change for humidity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDPRC', dim1name='column', &
      dim2name='month',long_name='accumulated change for precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDPRC,datrc_2d)    
  else
    !print*,'TDPRC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDPRC,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDPRC', dim1name='column', &
      dim2name='month',long_name='accumulated change for precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDIRI', dim1name='column', &
      dim2name='month',long_name='accumulated change for irrigation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDIRI,datrc_2d)    
  else
    !print*,'TDIRI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDIRI,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDIRI', dim1name='column', &
      dim2name='month',long_name='accumulated change for irrigation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCN4', dim1name='column', &
      dim2name='month',long_name='accumulated change for NH4 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDCN4,datrc_2d)    
  else
    !print*,'TDCN4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDCN4,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCN4', dim1name='column', &
      dim2name='month',long_name='accumulated change for NH4 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCNO', dim1name='column', &
      dim2name='month',long_name='accumulated change for NO3 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDCNO,datrc_2d)    
  else
    !print*,'TDCNO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDCNO,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDCNO', dim1name='column', &
      dim2name='month',long_name='accumulated change for NO3 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKC4', dim1name='column', &
      dim2name='levsoi',long_name='Ca-NH4 Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKC4,datrc_2d)    
  else
    !print*,'GKC4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKC4,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKC4', dim1name='column', &
      dim2name='levsoi',long_name='Ca-NH4 Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCH', dim1name='column', &
      dim2name='levsoi',long_name='Ca-H Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCH,datrc_2d)    
  else
    !print*,'GKCH'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCH,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCH', dim1name='column', &
      dim2name='levsoi',long_name='Ca-H Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCA', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Al Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCA,datrc_2d)    
  else
    !print*,'GKCA'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCA,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='GKCA', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Al Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCM', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Mg Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCM,datrc_2d)    
  else
    !print*,'GKCM'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCM,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCM', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Mg Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCN', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Na Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCN,datrc_2d)    
  else
    !print*,'GKCN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCN,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCN', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Na Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCK', dim1name='column', &
      dim2name='levsoi',long_name='Ca-K Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCK,datrc_2d)    
  else
    !print*,'GKCK'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCK,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCK', dim1name='column', &
      dim2name='levsoi',long_name='Ca-K Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  
  call destroy(datrc_3d)
  call destroy(datrc_4d)
  call destroy(datrc_5d)

  END SUBROUTINE restartnc_soil


  !-----------------------------------------------------------------------
  subroutine restFile_getfile( file, path )
    !
    ! !DESCRIPTION:
    ! Determine and obtain netcdf restart file
    !
    ! !USES:
    use EcoSIMConfig     , only : restartFileFullPath, nsrest, nsrStartup, nsrBranch 
    use EcoSIMConfig     , only : nsrContinue, case_name,  brnch_retain_casename
    use fileutil         , only : getfil
    !
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
    !
    ! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    
    ! Branch run: 
    ! Restart file pathname is obtained from namelist "restartFileFullPath"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)
    elseif (nsrest==nsrBranch) then
       length = len_trim(finidat)
       if (finidat(length-2:length) == '.nc') then
          path = trim(finidat) 
       else
          path = trim(finidat) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .elm makes this more robust
       ctest = 'xx.'//trim(case_name)//'.ecosim'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
!          if (masterproc) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(case_name), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
!          end if
          call endrun(msg=errMsg(__FILE__, __LINE__)) 
       end if
    end if

  end subroutine restFile_getfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

  use ncdio_pio
  implicit none
  character(len=*),  intent(in) :: flag ! flag to specify read or write
  character(len=*),  intent(in) :: file ! filename
  type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

!       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',etimer%get_nstep()
          write(iulog,*)
!       end if+
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

 !      if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
 !      end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

  !-----------------------------------------------------------------------
  subroutine restFile_closeRestart( file )
    !
    ! !DESCRIPTION:
    ! Close restart file and write restart pointer file if
    ! in write mode, otherwise just close restart file if in read mode
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS:
    character(len=*) , intent(in) :: file  ! local output filename
    !
    ! !CALLED FROM:
    ! subroutine restart in this module
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                   !index
    !-----------------------------------------------------------------------

!    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
!    end if

  end subroutine restFile_closeRestart

  !-----------------------------------------------------------------------
  subroutine restFile_write_pfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use EcoSIMConfig, only : rpntdir, rpntfil
    use fileutil , only : relavu
    use fileutil , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    character(len=datestrlen) :: curr_date
    !-----------------------------------------------------------------------

!    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )
       curr_date =etimer%get_calendar()
       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
!    end if

  end subroutine restFile_write_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_read_pfile( pnamer )
    !
    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    !
    ! !USES:
    use fileutil    , only : opnfil, getavu, relavu
    use EcoSIMConfig, only : rpntfil, rpntdir, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [restartFileFullPath] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.

!    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
!    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

!    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
!    end if

  end subroutine restFile_read_pfile  

  !-----------------------------------------------------------------------
  subroutine restFile_read( bounds, file)

  implicit none
  ! !ARGUMENTS:
  character(len=*)               , intent(in)    :: file  ! output netcdf restart file
  type(bounds_type)              , intent(in)    :: bounds  

! !LOCAL VARIABLES:
  type(file_desc_t) :: ncid         ! netcdf id
  integer           :: nc
  integer :: i

  ! Open file
  print*,trim(file)
  call restFile_open( flag='read', file=file, ncid=ncid )
  print*,'restFile_open Successfully'
  ! Read file

  call restFile_dimcheck( ncid )

  print*,'hist_restart_ncd'
  call hist_restart_ncd (bounds, ncid, flag='read')

  print*,'Do error checking on file'
  
  call restFile_check_consistency(bounds, ncid)
  
  call timemgr_restart_io( ncid, flag='read')
  
  call restartnc(ncid,flag='read')
  
  ! Close file 
  call restFile_close( ncid )
!    if (masterproc) then
      write(iulog,'(72a1)') ("-",i=1,60)
      write(iulog,*) 'Successfully read restart data for restart run'
      write(iulog,*)
!    end if

  end subroutine restFile_read  

  !-----------------------------------------------------------------------
  subroutine restFile_check_consistency(bounds, ncid)
  implicit none
  ! !ARGUMENTS:
  type(bounds_type), intent(in)    :: bounds  ! bounds
  type(file_desc_t), intent(inout) :: ncid    ! netcdf id
  character(len=*), parameter :: subname=trim(mod_filename)//'::restFile_check_consistency'

  print*,subname
  end subroutine restFile_check_consistency

  !-----------------------------------------------------------------------
  subroutine restFile_dimcheck( ncid )
    !
    ! !DESCRIPTION:
  implicit none 
  ! !ARGUMENTS:
  type(file_desc_t), intent(inout) :: ncid

  if ( .not. cold_run() )then
!      call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump, nCohorts=numCohort)
!      call check_dim(ncid, nameg, numg)
!      call check_dim(ncid, namet, numt)
!      call check_dim(ncid, namel, numl)
!      call check_dim(ncid, namec, numc)
!      call check_dim(ncid, namep, nump)
  end if

!  call check_dim(ncid, 'levsno'  , nlevsno)
!  call check_dim(ncid, 'levgrnd' , nlevgrnd)
!  call check_dim(ncid, 'levurb'  , nlevurb)
!  call check_dim(ncid, 'levlak'  , nlevlak) 

  call check_dim(ncid,'datestrlen',datestrlen)  
  call check_dim(ncid,'rootaxs',MaxNumRootAxes)
  call check_dim(ncid,'nodes',MaxNodesPerBranch)
  call check_dim(ncid,'nodes1',MaxNodesPerBranch+1)
  call check_dim(ncid,'cansecz',NumOfLeafZenithSectors)
  call check_dim(ncid,'fertN',trc_confs%NFertNitro)
  call check_dim(ncid,'fertNb',trc_confs%NFertNitrob)
  call check_dim(ncid,'automicb',NumMicrobAutrophCmplx)
  call check_dim(ncid,'nlbiomcp',nlbiomcp)
  call check_dim(ncid, 'levsoi', JZ)
  call check_dim(ncid, 'levsoi1', JZ+1)
  call check_dim(ncid, 'levsno',  JS)
  call check_dim(ncid, 'levcan',NumOfCanopyLayers)
  call check_dim(ncid, 'levcan1',NumOfCanopyLayers+1)
  call check_dim(ncid, 'npfts',  JP)
  call check_dim(ncid, 'nbranches',MaxNumBranches)
  call check_dim(ncid, 'ngrstages',NumGrowthStages)
  call check_dim(ncid, 'elmnts',NumPlantChemElms)
  call check_dim(ncid, 'nkinecmp',jsken)
  call check_dim(ncid, 'nomcomplx',jcplx)
  call check_dim(ncid, 'sdim',3)   !grid dimension
  call check_dim(ncid,'hetrmicb',NumHetetrMicCmplx)
  call check_dim(ncid,'rootyps'  , pltpar%jroots)
  call check_dim(ncid,'xtracers' , trc_confs%nxtracers)

  if(salt_model)then
    call check_dim(ncid,'satracers', trc_confs%NSaltTracers)
  endif  
  call check_dim(ncid,'ptracers' , trc_confs%NPrecipTracers)
  call check_dim(ncid,'NnutrientTracers', trc_confs%NnutrientTracers)
  call check_dim(ncid,'gastrcs'  , trc_confs%NGasTracers)
  call check_dim(ncid,'soltrcs'  , trc_confs%NSolutTracers)

  call check_dim(ncid , 'string_length', 64        )
  call check_dim(ncid , 'month'   , 12        )
  end subroutine restFile_dimcheck
  !-----------------------------------------------------------------------
  character(len=256) function restFile_filename( rdate )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
  implicit none
  !
  ! !ARGUMENTS:
  character(len=*), intent(in) :: rdate   ! input date for restart file name 
  !-----------------------------------------------------------------------

  restFile_filename = "./"//trim(case_name)//".ecosim"//trim(inst_suffix)//&
        ".r."//trim(rdate)//".nc"
  call strip_null(restFile_filename)      
  call strip_space(restFile_filename)

!    if (masterproc) then
    write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
!    end if

  end function restFile_filename  

  !-----------------------------------------------------------------------
  subroutine restFile_write( bounds, file,  rdate, noptr)

  implicit none
  type(bounds_type), intent(in) :: bounds
  character(len=*)               , intent(in)    :: file             ! output netcdf restart file  
  character(len=*)               , intent(in), optional :: rdate
  logical                        , intent(in), optional :: noptr     ! if should NOT write to the restart pointer file  
   ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! --------------------------------------------
    ! Open restart file
    ! --------------------------------------------

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset ( ncid )

    ! Define restart file variables

    call timemgr_restart_io(ncid, flag='define')
    
    if (present(rdate)) then       
      call hist_restart_ncd (bounds, ncid, flag='define', rdate=rdate )
    end if

    call restartnc(ncid,flag='define')

    call restFile_enddef( ncid )

    ! --------------------------------------------
    ! Write restart file variables
    ! --------------------------------------------
    
    call timemgr_restart_io( ncid, flag='write' )

    call hist_restart_ncd (bounds, ncid, flag='write' )

    call restartnc(ncid,flag='write')

    ! --------------------------------------------
    ! Close restart file and write restart pointer file
    ! --------------------------------------------
    
    call restFile_close( ncid )
    call restFile_closeRestart( file )

    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

!    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',etimer%get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
!    end if
    
  end subroutine restFile_write

  !=========================================================================================

  subroutine timemgr_restart_io( ncid, flag )

    !---------------------------------------------------------------------------------
    ! Read/Write information needed on restart to a netcdf file. 
    use ncdio_pio, only: ncd_int, var_desc_t, file_desc_t

    use EcosimConst, only : secspday
    !
    implicit none
    ! Arguments
    type(file_desc_t), intent(inout) :: ncid  ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
    !
    ! Local variables
    character(len=*), parameter :: sub = trim(mod_filename)//'::'//'timemgr_restart_io'
    integer :: rc                  ! return code
    logical :: readvar             ! determine if variable is on initial file
    character(len=datestrlen) :: curr_date   ! date of data in restart file, YYYYMMDDHHMMSS
    integer :: rst_caltype         ! calendar type
    integer :: isecspday
    integer :: varid
    integer, parameter :: noleap = 1
    integer, parameter :: gregorian = 2
    
    integer :: rst_start_ymd
    integer :: rst_ref_ymd
    integer :: rst_curr_ymd
    integer :: rst_start_tod
    integer :: rst_ref_tod
    integer :: rst_curr_tod
    integer :: rst_step_sec
    !---------------------------------------------------------------------------------

    if (flag == 'write') then
       rst_step_sec  = etimer%get_step_size()
       curr_date =etimer%get_calendar()
       read(start_date,'(I8)')rst_start_ymd       !the start time info of the simulation
       read(ref_date,'(I8)')rst_ref_ymd           !the referene time info of the simulation
       read(curr_date,'(I8)')rst_curr_ymd         !current date info 
       rst_start_tod=get_tod(start_date(9:14))
       rst_ref_tod  =get_tod(ref_date(9:14))
       rst_curr_tod =get_tod(curr_date(9:14))

       call ncd_putvar(ncid,'ref_date',ref_date)
       call ncd_putvar(ncid,'curr_date',curr_date)
       call ncd_putvar(ncid,'start_date',start_date)

    elseif(flag=='define')then

      call ncd_defvar(ncid,varname='ref_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='reference date of the simulation')
        
      call ncd_defvar(ncid,varname='start_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='starting date of the simulation')
        
      call ncd_defvar(ncid,varname='curr_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='current date of the simulation')

    elseif(flag=='read')then

    end if
    
    isecspday=secspday

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_step_sec', &
         long_name='seconds component of timestep size', units='sec',         &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                &
         interpinic_flag='skip',  data=rst_step_sec)
    
    if ((flag == 'read') .and. ( rst_step_sec < 0 .or. rst_step_sec > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_step_sec out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_ymd', &
         long_name='start date', units='YYYYMMDD', fill_value=uninit_int,     &
         interpinic_flag='skip',  data=rst_start_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_tod', &
         long_name='start time of day', units='sec',                           &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_start_tod)

    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_strart_tod out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_ymd',  &
         long_name='reference date', units='YYYYMMDD', fill_value=uninit_int, &
         interpinic_flag='skip',  data=rst_ref_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_tod',  &
         long_name='reference time of day', units='sec',                       &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_ref_tod)

    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_ref_tod out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_ymd', &
         long_name='current date', units='YYYYMMDD', fill_value=uninit_int,   &
         interpinic_flag='skip',  data=rst_curr_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_tod', &
         long_name='current time of day', units='sec',                         &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_curr_tod)

    if ((flag == 'read') .and. ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_ref_ymd out of range',__LINE__)
    end if

  end subroutine timemgr_restart_io
  !------------------------------------------------------------------------
  subroutine restFile_dimset( ncid )

  ! !ARGUMENTS:
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  !
  ! !LOCAL VARIABLES:
  integer :: dimid               ! netCDF dimension id
  integer :: numg                ! total number of gridcells across all processors
  integer :: numt                ! total number of topounits across all processors
  integer :: numc                ! total number of columns across all processors
  integer :: nump                ! total number of pfts across all processors

  integer :: ier                 ! error status
  integer :: strlen_dimid        ! string dimension id
  character(len=  8) :: curdate  ! current date
  character(len=  8) :: curtime  ! current time
  character(len=256) :: str
  character(len= 32), parameter :: subname=trim(mod_filename)//'::restFile_dimset' ! subroutine name
  !------------------------------------------------------------------------

  numg=bounds%ngrid;numt=bounds%ntopou;numc=bounds%ncols;nump=bounds%npfts

    ! Define dimensions

  call ncd_defdim(ncid , nameg      , numg           ,  dimid)
  call ncd_defdim(ncid , namet      , numt           ,  dimid)
  call ncd_defdim(ncid , namec      , numc           ,  dimid)
  call ncd_defdim(ncid , namep      , nump           ,  dimid)
  
  ! "level" dimensions
  call ncd_defdim(ncid,'datestrlen',datestrlen,dimid)  
  call ncd_defdim(ncid,'rootaxs',MaxNumRootAxes,dimid)
  call ncd_defdim(ncid,'nodes',MaxNodesPerBranch,dimid)
  call ncd_defdim(ncid,'nodes1',MaxNodesPerBranch+1,dimid)  
  call ncd_defdim(ncid,'cansecz',NumOfLeafZenithSectors,dimid)
  call ncd_defdim(ncid,'fertN',trc_confs%NFertNitro, dimid)
  call ncd_defdim(ncid,'fertNb',trc_confs%NFertNitrob, dimid)
  call ncd_defdim(ncid,'automicb',NumMicrobAutrophCmplx, dimid)
  call ncd_defdim(ncid,'nlbiomcp',nlbiomcp,dimid)
  call ncd_defdim(ncid,'nlivehetermicb',NumLiveHeterBioms,dimid)
  call ncd_defdim(ncid,'nliveautomicb',NumLiveAutoBioms,dimid)
  call ncd_defdim(ncid, 'levsoi', JZ, dimid)
  call ncd_defdim(ncid, 'levsoi1', JZ+1, dimid)
  call ncd_defdim(ncid, 'levsno',  JS,dimid)
  call ncd_defdim(ncid, 'levcan',NumOfCanopyLayers,dimid)
  call ncd_defdim(ncid, 'levcan1',NumOfCanopyLayers+1,dimid)
  call ncd_defdim(ncid, 'npfts',  JP,dimid)
  call ncd_defdim(ncid, 'nbranches',MaxNumBranches,dimid)
  call ncd_defdim(ncid, 'ngrstages',NumGrowthStages,dimid)
  call ncd_defdim(ncid, 'elmnts',NumPlantChemElms,dimid)
  call ncd_defdim(ncid, 'nkinecmp',jsken,dimid)
  call ncd_defdim(ncid, 'nomcomplx',jcplx,dimid)
  call ncd_defdim(ncid, 'nlitromcomplx',micpar%NumOfLitrCmplxs,dimid)
  call ncd_defdim(ncid, 'sdim',3,dimid)   !grid dimension
  call ncd_defdim(ncid,'hetrmicb',NumHetetrMicCmplx,dimid)
  call ncd_defdim(ncid,'rootyps'  , pltpar%jroots, dimid)
  call ncd_defdim(ncid,'xtracers' , trc_confs%nxtracers, dimid)
  call ncd_defdim(ncid,'ndoms',trc_confs%NDOMS,dimid)
  if(salt_model)call ncd_defdim(ncid,'satracers', trc_confs%NSaltTracers, dimid)
  call ncd_defdim(ncid,'ptracers' , trc_confs%NPrecipTracers, dimid)
  call ncd_defdim(ncid,'NnutrientTracers', trc_confs%NnutrientTracers, dimid)
  call ncd_defdim(ncid,'gastrcs'  , trc_confs%NGasTracers, dimid)
  call ncd_defdim(ncid,'soltrcs'  , trc_confs%NSolutTracers, dimid)

  call ncd_defdim(ncid , 'string_length', 64        ,  dimid)
  call ncd_defdim(ncid , 'month'   , 12        ,  dimid)


  if (do_budgets) then
  end if

    ! Define global attributes

!  call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
  call etimer%getdatetime(curdate, curtime)
  str = 'created on ' // curdate // ' ' // curtime
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(case_name)),subname)
!  call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
!  call ncd_putatt(ncid, NCD_GLOBAL, 'flanduse_timeseries', trim(get_flanduse_timeseries()))
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'title', 'ECOSIM Restart information'),subname)

!  call restFile_add_ipft_metadata(ncid)
!  call restFile_add_icol_metadata(ncid)

  end subroutine restFile_dimset  
  !-----------------------------------------------------------------------
  subroutine restFile_close( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
    !-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

!------------------------------------------------------------------------------------------

  integer function get_tod(hhmmss)
  !return the current time of day by seconds

  implicit none
  character(len=6) :: hhmmss
  integer :: hh,mm,ss

  read(hhmmss,'(I2,I2,I2)')hh,mm,ss
  get_tod=hh*3600+mm*60+ss

  end function get_tod

  !-----------------------------------------------------------------------
  subroutine restFile_enddef( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !-----------------------------------------------------------------------

    call ncd_enddef(ncid)

  end subroutine restFile_enddef

  !-----------------------------------------------------------------------
  subroutine get_restart_date(curr_date)
  !
  ! DESCRIPTION
  ! read the simulation date of restart file
  use EcoSIMConfig, only : rpntdir, rpntfil
  use fileutil , only : relavu
  use fileutil , only : getavu, opnfil
  implicit none
  character(len=datestrlen), intent(out) :: curr_date
  character(len=256) :: filename  ! local file name
  character(len=256) :: fnamer    ! restart file name
  integer :: nio
  integer :: ii,year,mon,day,tod,iostat
  nio = getavu()
  filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
  
  call opnfil( filename, nio, 'f' )  
  read(nio,'(A)',iostat=iostat) fnamer
  if (iostat == 0) then
    write(iulog,*) 'use restart file: ', fnamer
  else if (iostat > 0) then
    write(iulog,*) 'Error occurred during read operation: iostat = ', iostat
    call endrun(msg=errMsg(__FILE__, __LINE__)) 
  else
    write(iulog,*) 'End of file reached'
    call endrun(msg=errMsg(__FILE__, __LINE__)) 
  end if
  call relavu( nio )
  
  !get current date
  ii=2
  do
    ii=ii+1
    if(fnamer(ii-2:ii)=='.r.')exit
  enddo
  read(fnamer(ii+1:),'(I4)')year
  ii=ii+5
  read(fnamer(ii+1:),'(I2)')mon
  ii=ii+3
  read(fnamer(ii+1:),'(I2)')day
  ii=ii+3
  read(fnamer(ii+1:),'(I6)')tod
  write(curr_date,'(I4.4,I2.2,I2.2,I6.6)')year,mon,day,tod
  end subroutine get_restart_date
end module restartMod
