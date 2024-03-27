module RestartMod
!
! DESCRIPTION
! code to write restart/check point files
! code to read/write restart files
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use data_const_mod    , only : spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL  
  use EcoSIMConfig      , only : jcplx=> jcplxc, NumMicbFunGroups=> NumMicbFunGroups,nlbiomcp=>NumLiveMicrbCompts
  use EcoSIMConfig      , only : ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc, cold_run
  use EcoSIMConfig      , only : inst_suffix,ref_date,start_date, ctitle
  use EcoSIMConfig      , only : case_name,hostname,version,source,username
  use EcoSIMCtrlDataType, only : iYearCurrent
  use EcoSiMParDataMod  , only : micpar,pltpar
  use TracerIDMod       , only : trc_confs
  use EcoSIMCtrlMod     , only : etimer,do_budgets
  use restUtilMod  
  use abortutils        , only : endrun,destroy
  use HistFileMod       , only : hist_restart_ncd  
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
  use PlantMngmtDataType
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
  subroutine restartnc(ncid,flag)
  implicit none
  character(len=*), intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id

  call wouts(ncid,flag)

  CALL WOUTP(ncid,flag)
  end subroutine restartnc
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
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,tod
    filer = restFile_filename(rdate=rdate)
    call restFile_write(bounds, filer,rdate=rdate)
  endif
  end subroutine restFile

!------------------------------------------------------------------------------------------
  SUBROUTINE woutp(ncid,flag)  
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
  real(r8), target :: datrc_2d(bounds%begc:bounds%endc,1:JZ+1)
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
    datpr3 => datrp_3d(1:npfts,1:jsken,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadKCompElms_pft', dim1name='pft',dim2name='nkinecmp',&
       dim3name='elmnts',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadKCompElms_pft,datrp_3d)
  else
    !print*,'StandDeadKCompElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadKCompElms_pft,datrp_3d)
    datpr3 => datrp_3d(1:npfts,1:jsken,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='StandDeadKCompElms_pft', dim1name='pft',dim2name='nkinecmp',&
       dim3name='elmnts',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)

  endif

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlanting_pft', dim1name='pft',&
       long_name='year of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlanting_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)

  else
    !print*,'iYearPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iYearPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
    call restartvar(ncid, flag, varname='iYearPlantHarvest_pft', dim1name='pft',&
       long_name='year of harvest', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    

  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iDayPlantHarvest_pft', dim1name='pft',&
     long_name='day of harvest', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
  else 
    !print*,'iDayPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iDayPlantHarvest_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NGTopRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
  else
    !print*,'NGTopRootLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NGTopRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantShootState_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)

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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantRootState_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MainBranchNum_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitPlant_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call restartvar(ncid, flag, varname='NumRootAxes_pft', dim1name='pft',&
     long_name='root primary axis number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumRootAxes_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_pft,datip_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NIXBotRootLayer_rpft,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='PlantExudChemElmCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantExudChemElmCum_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PlantExudChemElmCum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantExudChemElmCum_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)  
    call restartvar(ncid, flag, varname='PlantExudChemElmCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='LitrfalStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element LitrFall', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LitrfalStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='PlantN2FixCum_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantN2FixCum_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PlantN2FixCum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PlantN2FixCum_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantN2FixCum_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossResp_pft', dim1name='pft',&
     long_name='total plant respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossResp_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='ETCanopy_pft', dim1name='pft',&
     long_name='total transpiration (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ETCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'ETCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ETCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ETCanopy_pft', dim1name='pft',&
     long_name='total transpiration  (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossCO2Fix_pft', dim1name='pft',&
     long_name='total gross CO2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrossCO2Fix_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TCelciusCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='TCG', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TCG,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TCG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TCG,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCG', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKG', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TKG,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'TKG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TKG,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKG', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           

  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTgrowCanP', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    call cppft(flag,NHW,NHE,NVN,NVS,NP,fTgrowCanP,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'fTgrowCanP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,fTgrowCanP,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTgrowCanP', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStalkC_pft', dim1name='pft',&
     long_name='canopy active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStalkC_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyWater_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyOsmo_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSICanopyTurg_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='RootBiomCPerPlant_pft', dim1name='pft',&
     long_name='root C per plant', units='g g-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootBiomCPerPlant_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  

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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafArea_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmnt_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'EcoHavstElmnt_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmnt_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmntCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
     call cppft(flag,NHW,NHE,NVN,NVS,NP,EcoHavstElmntCum_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='CanopyPlusNoduRespC_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyPlusNoduRespC_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'CanopyPlusNoduRespC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyPlusNoduRespC_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyPlusNoduRespC_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NetCumElmntFlx2Plant_pft', dim1name='pft',dim2name='elmnts',&
     long_name='effect of canopy element status on seed set', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NetCumElmntFlx2Plant_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='NH3EmiCum_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3EmiCum_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'NH3EmiCum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3EmiCum_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3EmiCum_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  

  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElmsCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SurfLitrfalStrutElmsCum_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'SurfLitrfalStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,SurfLitrfalStrutElmsCum_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElmsCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI', dim1name='pft',&
     long_name='initial plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PPI,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PPI'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PPI,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI', dim1name='pft',&
     long_name='initial plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX', dim1name='pft',&
     long_name='plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PPX,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'PPX'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PPX,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX', dim1name='pft',&
     long_name='plant population', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='standing dead element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StandDeadStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyHeight_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
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
    call restartvar(ncid, flag, varname='WatByPCanopy', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            
    call cppft(flag,NHW,NHE,NVN,NVS,NP,WatByPCanopy,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'WatByPCanopy'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,WatByPCanopy,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='WatByPCanopy', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  
  
  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ClumpFactor,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)  
  else
    !print*,'ClumpFactor'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ClumpFactor,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
  else
    !print*,'CO2ByFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CO2ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CH4ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
  else
    !print*,'CH4ByFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CH4ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)      
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP,O2ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'O2ByFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,O2ByFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3byFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'NH3byFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,NH3byFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,N2ObyFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'N2ObyFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,N2ObyFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PO4byFire_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'PO4byFire_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PO4byFire_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='MatureGroup_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,MatureGroup_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
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
    call restartvar(ncid, flag, varname='ShootNodeNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootNodeNumber_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
  else
    !print*,'ShootNodeNumber_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootNodeNumber_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)        
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNum2InitFloral_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNum2InitFloral_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumberAtAnthesis_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)      
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NumOfLeaves_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)        
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafNumberAtFloralInit_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)        
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)        
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'ReprodNodeNumNormByMatrgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,TotalNodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'TotalNodeNumNormByMatgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,TotalNodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=NumActivePlants,&
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4Leafout_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LeafOff_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LenthenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4LenthenPhotoPeriod_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LenthenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,&
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4ShortenPhotoPeriod_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
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
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitLeafOut_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doInitLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doInitLeafOut_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeafOut_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doPlantLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeafOut_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'doPlantLeaveOff_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Prep4Literfall_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Prep4Literfall_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Prep4Literfall_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)            
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'Hours4LiterfalAftMature_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantBranchState_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantBranchState_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'iPlantBranchState_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantBranchState_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantBranchState_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'BranchNumber_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,BranchNumber_brch,datip_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KLeafNumber_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)          
  else
    !print*,'KLEAF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KLeafNumber_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)    
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'KVSTG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'KLowestGroLeafNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    dat3pr => datip_3d(1:npfts,1:NumGrowthStages,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantCalendar_brch', dim1name='pft',dim2name='ngrstages', &
     dim3name='nbranches',long_name='plant growth stage', units='none', &
     interpinic_flag='skip', data=dat3pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,iPlantCalendar_brch,datip_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='VHeatCapCanP', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,VHeatCapCanP,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'VHeatCapCanP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,VHeatCapCanP,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='VHeatCapCanP', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DTKC', dim1name='pft',&
     long_name='change in canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,DTKC,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'DTKC'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,DTKC,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DTKC', dim1name='pft',&
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
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'LeafStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LeafStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetioleStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
  else
    !print*,'PetioleStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)     
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetioleStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,NodulStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedDepth_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
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
    call restartvar(ncid, flag, varname='HoursCanopyPSITooLow', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursCanopyPSITooLow,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
  else
    !print*,'HoursCanopyPSITooLow'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,HoursCanopyPSITooLow,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HoursCanopyPSITooLow', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CHILL', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CHILL,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
  else
    !print*,'CHILL'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CHILL,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CHILL', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracRadPARbyCanopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,FracRadPARbyCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)     
  else
    !print*,'FracRadPARbyCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,FracRadPARbyCanopy_pft,datrp_1d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)       
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracRadPARbyCanopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)   

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootStrutElms_pft,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)       
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkBiomassC_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulNonstElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ShootStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyNodulStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetoleStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,StalkRsrvElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,HuskStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,EarStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PotentialSeedSites_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PotentialSeedSites_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PotentialSeedSites_brch,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SeedNumSet_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch grain number', units='# d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SeedNumSet_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,GrainSeedBiomCMean_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaDying_brch,datrp_2d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemob_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleChemElmRemobFlx_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,SenecStalkStrutElms_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='CanopyLeafALyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafALyr_pft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyLeafALyr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafALyr_pft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafALyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemArea_lpft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_lpft,datrp_2d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyStemArea_lpft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_lpft,datrp_2d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr2 => datrp_2d(1:npfts,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemArea_lpft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL4_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C4 content in C4 photosynthesi', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL4_node,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CPOOL3_node,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassCO2BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CMassHCO3BundleSheath_node,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaNode_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'ARLF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafAreaNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinCNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafProteinCNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinCNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleLengthNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleLengthNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PetioleLengthNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleLengthNode_brch,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleLengthNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif 

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'WGSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleElmntNode_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)      
  endif  
 
  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleProteinCNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'WSSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PetioleProteinCNode_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetioleProteinCNode_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightLive_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightLive_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'InternodeHeightLive_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightLive_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightLive_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='living internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightDying_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightDying_brch,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'InternodeHeightDying_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeHeightDying_brch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeHeightDying_brch', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeStrutElms_brch,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'InternodeStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,InternodeStrutElms_brch,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='InternodeStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafAreaByLayer_pft', dim1name='pft',dim2name='levcan',&
     dim3name='nodes',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafAreaByLayer_pft,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyLeafAreaByLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyLeafAreaByLayer_pft,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafAreaByLayer_pft', dim1name='pft',dim2name='levcan',&
     dim3name='nodes',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafChemElmByLayerNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='layer leaf element',&
     units='g d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafChemElmByLayerNode_brch,datrp_5d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'LeafChemElmByLayerNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,LeafChemElmByLayerNode_brch,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumOfCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafChemElmByLayerNode_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='layer leaf element',&
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
    call restartvar(ncid, flag, varname='CanopyStemArea_lbrch', dim1name='pft',dim2name='levcan',&
     dim3name='nbranches',long_name='plant canopy layer branch stem layer area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_lbrch,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'CanopyStemArea_lbrch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,CanopyStemArea_lbrch,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:NumOfCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyStemArea_lbrch', dim1name='pft',dim2name='levcan',&
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

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRootOSMO_vr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root osmotic water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootOSMO_vr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,PSIRootTurg_vr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stXNumL_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootLenDensPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootPoreVol_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootVH2O_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stRadius_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='Radius2ndRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Radius2ndRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'Radius2ndRoot_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Radius2ndRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='Radius2ndRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RootAreaPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root layer area per plant', units='m p-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootAreaPerPlant_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='AveLen2ndRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,AveLen2ndRoot_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'AveLen2ndRoot_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,AveLen2ndRoot_pvr,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='AveLen2ndRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RootRespPotent_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root respiration unconstrained by O2', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootRespPotent_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='RCO2A_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RCO2A_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RCO2A_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RCO2A_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RCO2A_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RCO2N_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RCO2N_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RCO2N_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RCO2N_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RCO2N_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='trcg_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
      dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP,trcg_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,trcs_rootml_pvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMycoActiveBiomC_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname=' PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP, PopuRootMycoC_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,' PopuRootMycoC_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP, PopuRootMycoC_pvr,datrp_3d,&
      NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname=' PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RootProteinC_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root layer protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootProteinC_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='ROXYP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,ROXYP,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'ROXYP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,ROXYP,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='ROXYP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNHP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNHP,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUNNHP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNHP,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNHP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNOP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNOP,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUNNOP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNOP,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNOP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP2P', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP2P,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUPP2P'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP2P,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP2P', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP1P', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP1P,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUPP1P'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP1P,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP1P', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNBP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNBP,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUNNBP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNBP,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNBP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNXP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNXP,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUNNXP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUNNXP,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUNNXP', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP2B', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP2B,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUPP2B'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP2B,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP2B', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP1B', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP1B,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RUPP1B'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RUPP1B,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RUPP1B', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RAutoRootO2Limter_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RAutoRootO2Limter_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'WFR'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RAutoRootO2Limter_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft) 
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='RAutoRootO2Limter_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname=' RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
     dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP, RootMycoNonstElms_rpvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,' RootMycoNonstElms_rpvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP, RootMycoNonstElms_rpvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname=' RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
     dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts, 1:pltpar%jroots,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root1stDepz_pft', dim1name='pft',dim2name='rootyps',&
     dim3name='levcan',long_name='primary root layer depth', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stDepz_pft,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    call restartvar(ncid, flag, varname='Root1stElm_raxs', dim1name='pft',dim2name='elmnts',&
     dim3name='rootyps',dim4name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stElm_raxs,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'Root1stElm_raxs'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stElm_raxs,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='Root1stElm_raxs', dim1name='pft',dim2name='elmnts',&
     dim3name='rootyps',dim4name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes number', units=' d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RTN2'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes number', units=' d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root1stLen_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer length primary axes', units='m d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root1stLen_rpvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
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
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer element primary axes', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RootMyco1stStrutElms_rpvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco1stStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer primary axes element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RootMyco2ndStrutElms_rpvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer secondary axes element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root2ndLen_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer length secondary axes', units='m d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndLen_pvr,datrp_4d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'Root2ndLen_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,Root2ndLen_pvr,datrp_4d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:NumOfCanopyLayers)
    call restartvar(ncid, flag, varname='Root2ndLen_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',dim4name='levcan',long_name='root layer length secondary axes', units='m d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
    call restartvar(ncid, flag, varname='RootNodulNonstElms_pvr', dim1name='pft',dim2name='elmnts',&
     dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulNonstElms_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RootNodulNonstElms_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulNonstElms_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)  
    datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
    call restartvar(ncid, flag, varname='RootNodulNonstElms_pvr', dim1name='pft',dim2name='elmnts',&
     dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
    call restartvar(ncid, flag, varname='RootNodulStrutElms_pvr', dim1name='pft',dim2name='elmnts',&
     dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulStrutElms_pvr,datrp_3d,NumActivePlants=NumActivePlants,IsPlantActive_pft=IsPlantActive_pft) 
  else
    !print*,'RootNodulStrutElms_pvr'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP,RootNodulStrutElms_pvr,datrp_3d,NumActivePlants=NumActivePlants,&
      IsPlantActive_pft=IsPlantActive_pft)   
    datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
    call restartvar(ncid, flag, varname='RootNodulStrutElms_pvr', dim1name='pft',dim2name='elmnts',&
     dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  call destroy(datrp_3d)
  call destroy(datrp_4d)  
  call destroy(datrp_5d)  
  call destroy(datip_3d)

  END SUBROUTINE woutp

!------------------------------------------------------------------------------------------

  SUBROUTINE wouts(ncid,flag)
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
          NumMicrobAutotrophCmplx,jcplx)
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

  call restartvar(ncid, flag, varname='TSEDOU', &
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=TSEDOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HEATIN', &
       long_name='total surface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HEATIN,missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='OXYGIN', &
       long_name='total surface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=OXYGIN, missing_value=spval, &
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

  call restartvar(ncid, flag, varname='ZN2GIN', &
       long_name='total surface N2 flux', units='g d-2', &
       interpinic_flag='skip', data=ZN2GIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CO2GIN', &
       long_name='total surface CO2 flux', units='g d-2', &
       interpinic_flag='skip', data=CO2GIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='VOLWOU', &
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=VOLWOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CEVAP', &
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=CEVAP, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CRUN', &
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=CRUN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HEATOU', &
       long_name='total subsurface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HEATOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='OXYGOU', &
       long_name='total subsurface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=OXYGOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TCOU', &
       long_name='total subsurface C flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TZOU', &
       long_name='total subsurface N flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU(ielmn), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TPOU', &
       long_name='total subsurface P flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU(ielmp), missing_value=spval, &
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
       interpinic_flag='skip', data=XESN(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='XZSN', &
       long_name='total LitrFall N', units='g d-2', &
       interpinic_flag='skip', data=XESN(ielmn), missing_value=spval, &
       fill_value=spval)

  
  call restartvar(ncid, flag, varname='XPSN', &
       long_name='total LitrFall P', units='g d-2', &
       interpinic_flag='skip', data=XESN(ielmp), missing_value=spval, &
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
    call restartvar(ncid, flag, varname='MaxCanopyHeight_grd', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,MaxCanopyHeight_grd,datrc_1d)
  else
    !print*,'MaxCanopyHeight_grd'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,MaxCanopyHeight_grd,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='MaxCanopyHeight_grd', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='PrecIntcptByCanG', dim1name='column',&
       long_name='net water transfer to whole grid canopy', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,PrecIntcptByCanG,datrc_1d)
  else
    !print*,'PrecIntcptByCanG'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PrecIntcptByCanG,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='PrecIntcptByCanG', dim1name='column',&
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
    call restartvar(ncid, flag, varname='SoiSurfRoughnesst0', dim1name='column',&
       long_name='initial soil surface roughness height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughnesst0,datrc_1d)
  else
    !print*,'SoiSurfRoughnesst0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughnesst0,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SoiSurfRoughnesst0', dim1name='column',&
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
    call restartvar(ncid, flag, varname='LWRadBySurf', dim1name='column',&
       long_name='longwave radiation emitted from ground surface (including snow and litter)', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf,datrc_1d)
  else
    !print*,'LWRadBySurf'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='LWRadBySurf', dim1name='column',&
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
    call restartvar(ncid, flag, varname='CanH2OHeldVg', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CanH2OHeldVg,datrc_1d) 
  else
    !print*,'CanH2OHeldVg'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanH2OHeldVg,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanH2OHeldVg', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='VWatStoreCapSurf', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf,datrc_1d) 
  else
    !print*,'VWatStoreCapSurf'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='VWatStoreCapSurf', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='URAIN', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,URAIN,datrc_1d) 
  else
    !print*,'URAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,URAIN,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='URAIN', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyLeafArea_grd', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_grd,datrc_1d) 
  else
    !print*,'CanopyLeafArea_grd'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_grd,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyLeafArea_grd', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='StemArea_grd', dim1name='column',&
       long_name='grid canopy stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_grd,datrc_1d) 
  else
    !print*,'StemArea_grd'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_grd,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='StemArea_grd', dim1name='column',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,PPT,datrc_1d) 
  else
    !print*,'PPT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PPT,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='PPT', dim1name='column',&
       long_name='total plant population', units='# d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='MaxVLWatByLitR', dim1name='column',&
       long_name='soil surface water retention capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,MaxVLWatByLitR,datrc_1d) 
  else
    !print*,'MaxVLWatByLitR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,MaxVLWatByLitR,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='MaxVLWatByLitR', dim1name='column',&
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
    call restartvar(ncid, flag, varname='SurfGasFlx', dim1name='column',&
       dim2name='gastrcs',long_name='total soil gas flux', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SurfGasFlx,datrc_2d) 
  else
    !print*,'UCH4G'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SurfGasFlx,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)          
    call restartvar(ncid, flag, varname='SurfGasFlx', dim1name='column',&
       dim2name='gastrcs',long_name='total soil gas flux', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CO2byFire_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_col,datrc_1d) 
  else
    !print*,'CO2byFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CO2byFire_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    

  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CH4byFire_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_col,datrc_1d) 
  else
    !print*,'CH4byFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CH4byFire_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='O2byFire_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_col,datrc_1d) 
  else
    !print*,'O2byFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='O2byFire_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='N2ObyFire_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_col,datrc_1d) 
  else
    !print*,'N2ObyFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_col,datrc_1d)
    datpr1 => datrc_1d             
    call restartvar(ncid, flag, varname='N2ObyFire_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='NH3byFire_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_col,datrc_1d) 
  else
    !print*,'NH3byFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='NH3byFire_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='PO4byFire_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_col,datrc_1d) 
  else
    !print*,'PO4byFire_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='PO4byFire_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='AmendCFlx_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,AmendCFlx_col,datrc_1d) 
  else
    !print*,'AmendCFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AmendCFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='AmendCFlx_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FertNFlx_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FertNFlx_col,datrc_1d) 
  else
    !print*,'FertNFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertNFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FertNFlx_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FerPFlx_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FerPFlx_col,datrc_1d) 
  else
    !print*,'FerPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FerPFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FerPFlx_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='UVOLO', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,UVOLO,datrc_1d) 
  else
    !print*,'UVOLO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,UVOLO,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='UVOLO', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='UEVAP', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,UEVAP,datrc_1d) 
  else
    !print*,'UEVAP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,UEVAP,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='UEVAP', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='URUN', dim1name='column',&
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,URUN,datrc_1d) 
  else
    !print*,'URUN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,URUN,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='URUN', dim1name='column',&
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
    call restartvar(ncid, flag, varname='HydroSufDONFlx_col', dim1name='column',&
       long_name='total surface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDONFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDONFlx_col', dim1name='column',&
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
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_col', dim1name='column',&
       long_name='total surface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDOPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_col', dim1name='column',&
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
    call restartvar(ncid, flag, varname='HydroSufDINFlx_col', dim1name='column',&
       long_name='total surface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDINFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDINFlx_col', dim1name='column',&
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
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_col', dim1name='column',&
       long_name='total surface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDIPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_col', dim1name='column',&
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
       dim2name='elements',long_name='total LitrFall C', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col,datrc_2d) 
  else
    !print*,'LiterfalOrgC_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col,datrc_2d)   
    call restartvar(ncid, flag, varname='LiterfalOrgM_col', dim1name='column',&
       dim2name='elements',long_name='total LitrFall C', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='UCOP', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,UCOP,datrc_1d) 
  else
    !print*,'UCOP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,UCOP,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='UCOP', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='UDRAIN', dim1name='column',&
       long_name='total water drainage below root zone', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,UDRAIN,datrc_1d) 
  else
    !print*,'UDRAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,UDRAIN,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='UDRAIN', dim1name='column',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DPNH4,datrc_1d) 
  else
    !print*,'DPNH4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPNH4,datrc_1d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DPNO3,datrc_1d) 
  else
    !print*,'DPNO3'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPNO3,datrc_1d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DPPO4,datrc_1d) 
  else
    !print*,'DPPO4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPPO4,datrc_1d)   
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
    call restartvar(ncid, flag, varname='USEDOU', dim1name='column',&
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,USEDOU,datrc_1d) 
  else
    !print*,'USEDOU'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,USEDOU,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='USEDOU', dim1name='column',&
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
    call restartvar(ncid, flag, varname='Eco_NBP_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_col,datrc_1d) 
  else
    !print*,'Eco_NBP_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NBP_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='VLitR', dim1name='column',&
       long_name='surface litter volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLitR,datrc_1d) 
  else
    !print*,'VLitR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLitR,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='VLitR', dim1name='column',&
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
    call restartvar(ncid, flag, varname='Eco_GPP_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_col,datrc_1d) 
  else
    !print*,'Eco_GPP_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_GPP_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_AutoR_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_col,datrc_1d) 
  else
    !print*,'Eco_AutoR_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_AutoR_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_NPP_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_col,datrc_1d) 
  else
    !print*,'Eco_NPP_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NPP_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_HR_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_col,datrc_1d) 
  else
    !print*,'Eco_HR_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_HR_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Canopy_Heat_Latent_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Canopy_Heat_Latent_col,datrc_1d) 
  else
    !print*,'Canopy_Heat_Latent_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Canopy_Heat_Latent_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Canopy_Heat_Latent_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Canopy_Heat_Sens_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Canopy_Heat_Sens_col,datrc_1d) 
  else
    !print*,'Canopy_Heat_Sens_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Canopy_Heat_Sens_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Canopy_Heat_Sens_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DayLenthCurrent', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthCurrent,datrc_1d) 
  else
    !print*,'DayLenthCurrent'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthCurrent,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DayLenthCurrent', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='DayLenthPrev', dim1name='column',&
       long_name='daylength of previous day', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev,datrc_1d) 
  else
    !print*,'DayLenthPrev'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev,datrc_1d)   
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='DayLenthPrev', dim1name='column',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable,datrc_1d) 
  else
    !print*,'ExtWaterTable'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable,datrc_1d)   
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
    datpr2 => datrc_2d(1:ncols,1:jcplx)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nomcomplx',&
       long_name='surface mic OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RC0,datrc_2d) 
  else
    !print*,'RC0',size(RC0,1),size(RC0,2)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RC0,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:jcplx)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nomcomplx',&
       long_name='surface mic OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='RC0ff', dim1name='column',&
       long_name='surface mic OM in autotroph complex', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RC0ff,datrc_1d) 
  else
    !print*,'RC0ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RC0ff,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='RC0ff', dim1name='column',&
       long_name='surface mic OM in autotroph complex', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE,datrc_1d) 
  else
    !print*,'VcumDrySnoWE'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow,datrc_1d) 
  else
    !print*,'VcumWatSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow,datrc_1d) 
  else
    !print*,'VcumIceSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI', dim1name='column',&
       long_name='total snowpack volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI,datrc_1d) 
  else
    !print*,'VcumSnoDWI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI', dim1name='column',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth,datrc_1d) 
  else
    !print*,'SnowDepth'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth,datrc_1d) 
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE,datrc_2d) 
  else
    !print*,'VLDrySnoWE'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE,datrc_2d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow,datrc_2d) 
  else
    !print*,'VLIceSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow,datrc_2d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow,datrc_2d) 
  else
    !print*,'VLWatSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLWatSnow', dim1name='column',dim2name='levsno',&
       long_name='snow water volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWI', dim1name='column',dim2name='levsno',&
       long_name='total snow volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWI,datrc_2d) 
  else
    !print*,'VLSnoDWI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWI,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWI', dim1name='column',dim2name='levsno',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,SnoDensL,datrc_2d) 
  else
    !print*,'SnoDensL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnoDensL,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnoDensL', dim1name='column',dim2name='levsno',&
       long_name='snowpack density', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowLayerThick', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowLayerThick,datrc_2d) 
  else
    !print*,'SnowLayerThick'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowLayerThick,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowLayerThick', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow', dim1name='column',dim2name='levsno',&
       long_name='snowpack heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow,datrc_2d) 
  else
    !print*,'VLHeatCapSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow', dim1name='column',dim2name='levsno',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow,datrc_2d) 
  else
    !print*,'TKSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow,datrc_2d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow,datrc_2d) 
  else
    !print*,'TCSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow,datrc_2d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml,datrc_3d) 
  else
    !print*,'trcg_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml,datrc_3d)   
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
    call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml,datrc_3d) 
  else
    !print*,'trcn_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NnutrientTracers,1:JS)            
    call restartvar(ncid, flag, varname='trcn_solsml', dim1name='column',dim2name='NnutrientTracers',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SoilFracAsMacP', dim1name='column',dim2name='levsoi',&
       long_name='macropore fraction', units='m3/m3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP,datrc_2d) 
  else
    !print*,'SoilFracAsMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoilFracAsMacP', dim1name='column',dim2name='levsoi',&
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
    call restartvar(ncid, flag, varname='CumDepth2LayerBottom', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumDepth2LayerBottom,datrc_2d) 
  else
    !print*,'CumDepth2LayerBottom'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumDepth2LayerBottom,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumDepth2LayerBottom', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness,datrc_2d) 
  else
    !print*,'CumSoilThickness'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoiBulkDensityt0', dim1name='column',dim2name='levsoi',&
       long_name='initial bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0,datrc_2d) 
  else
    !print*,'SoiBulkDensityt0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='SoiBulkDensityt0', dim1name='column',dim2name='levsoi',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensity,datrc_2d) 
  else
    !print*,'SoiBulkDensity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensity,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoiBulkDensity', dim1name='column',dim2name='levsoi1',&
       long_name='soil bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='CORGC', dim1name='column',dim2name='levsoi1',&
       long_name='soil organic C content', units='g kg-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CORGC,datrc_2d) 
  else
    !print*,'CORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CORGC,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='CORGC', dim1name='column',dim2name='levsoi1',&
       long_name='soil organic C content', units='g kg-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    

  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,POROS,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,POROS,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity,datrc_2d) 
  else
    !print*,'FieldCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint,datrc_2d)   
  else
    !print*,'WiltPoint'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert,datrc_2d)   
  else
    !print*,'SatHydroCondVert'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn', dim1name='column',dim2name='levsoi',&
       long_name='soil horizontal saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn,datrc_2d)   
  else
    !print*,'SatHydroCondHrzn'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn', dim1name='column',dim2name='levsoi',&
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
    call restartvar(ncid, flag, varname='VLWatMicP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP,datrc_2d)   
  else
    !print*,'VLWatMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicP', dim1name='column',dim2name='levsoi1',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX,datrc_2d)     
  else
    !print*,'VLWatMicPX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX,datrc_2d)       
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP,datrc_2d)     
  else
    !print*,'VLiceMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP,datrc_2d)       
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP,datrc_2d)     
  else
    !print*,'VLsoiAirP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLsoiAirP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore air content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP,datrc_2d)     
  else
    !print*,'VLMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP,datrc_2d)     
  else
    !print*,'VLSoilMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)    
    call restartvar(ncid, flag, varname='trcx_solml', dim1name='column',dim2name='xtracers',&
       dim3name='levsoi1',long_name='exchangeable tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml,datrc_3d)     
  else
    !print*,'trcx_solml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)        
    call restartvar(ncid, flag, varname='trcx_solml', dim1name='column',dim2name='xtracers',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,TKS,datrc_2d)     
  else
    !print*,'TKS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKS,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TKS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity,datrc_2d)     
  else
    !print*,'VHeatCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacitySoilM', dim1name='column',dim2name='levsoi1',&
       long_name='soil solid heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacitySoilM,datrc_2d)     
  else
    !print*,'VHeatCapacitySoilM'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacitySoilM,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacitySoilM', dim1name='column',dim2name='levsoi1',&
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
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trc_gasml_vr,datrc_3d)       
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
    call restartvar(ncid, flag, varname='trc_soHml', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trc_soHml,datrc_3d)     
  else
    !print*,'trc_soHml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trc_soHml,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trc_soHml', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   

  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROXYF', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROXYF,datrc_2d)     
  else
    !print*,'ROXYF'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROXYF,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYF', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCO2F', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CO2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RCO2F,datrc_2d)     
  else
    !print*,'RCO2F'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCO2F,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCO2F', dim1name='column',dim2name='levsoi1',&
       long_name='net gaseous CO2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYL', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROXYL,datrc_2d)       
  else
    !print*,'ROXYL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROXYL,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYL', dim1name='column',dim2name='levsoi1',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,RCH4L,datrc_2d)       
  else
    !print*,'RCH4L'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4L,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4L', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROXYX,datrc_2d)      
  else
    !print*,'ROXYX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROXYX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNH4X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4X,datrc_2d)      
  else
    !print*,'RNH4X'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4X,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNH4X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO3X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3X,datrc_2d)      
  else
    !print*,'RNO3X'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3X,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO3X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2X,datrc_2d)      
  else
    !print*,'RNO2X'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,RNO2X,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2OX,datrc_2d)      
  else
    !print*,'RN2OX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2OX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RPO4X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RPO4X,datrc_2d)      
  else
    !print*,'RPO4X'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RPO4X,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RPO4X', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RP14X', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RP14X,datrc_2d)      
  else
    !print*,'RP14X'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RP14X,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RP14X', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNHBX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNHBX,datrc_2d)      
  else
    !print*,'RNHBX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNHBX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNHBX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN3BX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN3BX,datrc_2d)      
  else
    !print*,'RN3BX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN3BX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN3BX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2BX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2BX,datrc_2d)      
  else
    !print*,'RN2BX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2BX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2BX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RPOBX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RPOBX,datrc_2d)      
  else
    !print*,'RPOBX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RPOBX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RPOBX', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RP1BX', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RP1BX,datrc_2d)      
  else
    !print*,'RP1BX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RP1BX,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RP1BX', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)        
    call restartvar(ncid, flag, varname='ROQCX', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROQCX,datrc_3d)      
  else
    !print*,'ROQCX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROQCX,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='ROQCX', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='ROQAX', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial acetate uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROQAX,datrc_3d)      
  else
    !print*,'ROQAX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROQAX,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='ROQAX', dim1name='column',dim2name='nomcomplx',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP,datrc_2d)      
  else
    !print*,'VLWatMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP,datrc_2d)        
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
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP,datrc_2d)      
  else
    !print*,'VLiceMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLiceMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                
    call restartvar(ncid, flag, varname='VLMacP', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMacP,datrc_2d)      
  else
    !print*,'VLMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMacP,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='VLMacP', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='RTDNT', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,RTDNT,datrc_2d)      
  else
    !print*,'RTDNT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RTDNT,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='RTDNT', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                    
    call restartvar(ncid, flag, varname='CanopyHeightz_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightz_col,datrc_2d)      
  else
    !print*,'ZL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightz_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='CanopyHeightz_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)                        
    call restartvar(ncid, flag, varname='CanopyLAgrid_lyr', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLAgrid_lyr,datrc_2d)      
  else
    !print*,'CanopyLAgrid_lyr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLAgrid_lyr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyLAgrid_lyr', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemA_lyr', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemA_lyr,datrc_2d)      
  else
    !print*,'CanopyStemA_lyr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemA_lyr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemA_lyr', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='WGLFT', dim1name='column',dim2name='levcan',&
       long_name='total leaf mass', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,WGLFT,datrc_2d)      
  else
    !print*,'WGLFT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WGLFT,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)    
    call restartvar(ncid, flag, varname='WGLFT', dim1name='column',dim2name='levcan',&
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
    call restartvar(ncid, flag, varname='trcp_salml', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,trcp_salml,datrc_3d)      
  else
    !print*,'trcp_salml'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trcp_salml,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NPrecipTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trcp_salml', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)        
    call restartvar(ncid, flag, varname='CION', dim1name='column',dim2name='levcan',&
       long_name='solution ion concentratiom', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CION,datrc_2d)      
  else
    !print*,'CION'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CION,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumOfCanopyLayers)           
    call restartvar(ncid, flag, varname='CION', dim1name='column',dim2name='levcan',&
       long_name='solution ion concentratiom', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  


  IF(salt_model)THEN

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)           
      call restartvar(ncid, flag, varname='trcSalt_solml', dim1name='column',dim2name='satracers',&
        dim3name='levsoi1',long_name='soil aqueous salt content micropre', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml,datrc_3d)      
    else
      !print*,'trcSalt_solml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)               
      call restartvar(ncid, flag, varname='trcSalt_solml', dim1name='column',dim2name='satracers',&
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
        call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solsml,datrc_3d)      
    else
      !print*,'trcs_solsml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solsml,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JS)    
      call restartvar(ncid, flag, varname='trcs_solsml', dim1name='column',dim2name='satracers',&
         dim3name='levsno',long_name='snowpack salt dissolved tracers', units='mol d-2', &
         interpinic_flag='skip', data=datpr3, missing_value=spval, &
         fill_value=spval)  

    endif  

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml,datrc_3d)      
    else
      !print*,'trcSalt_soHml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml,datrc_3d)          
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
    endif  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROXYS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROXYS,datrc_4d)      
  else
    !print*,'ROXYS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROXYS,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)        
    call restartvar(ncid, flag, varname='ROXYS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)        
    call restartvar(ncid, flag, varname='RVMX4', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NH4 uptake non-band unconstrained by NH4',&
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX4,datrc_4d)      
  else
    !print*,'RVMX4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX4,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)            
    call restartvar(ncid, flag, varname='RVMX4', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NH4 uptake non-band unconstrained by NH4',&
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX3', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX3,datrc_4d)      
  else
    !print*,'RVMX3'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX3,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX3', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX2', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX2,datrc_4d)      
  else
    !print*,'RVMX2'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX2,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX2', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX1', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX1,datrc_4d)      
  else
    !print*,'RVMX1'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX1,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX1', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB4', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB4,datrc_4d)      
  else
    !print*,'RVMB4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB4,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB4', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB3', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB3,datrc_4d)      
  else
    !print*,'RVMB3'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB3,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB3', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB2', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB2,datrc_4d)      
  else
    !print*,'RVMB2'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB2,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB2', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHO,datrc_4d)      
  else
    !print*,'RINHO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHO,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOO,datrc_4d)      
  else
    !print*,'RINOO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOO,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPOO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPOO,datrc_4d)      
  else
    !print*,'RIPOO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPOO,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPOO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHB', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHB,datrc_4d)      
  else
    !print*,'RINHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHB,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHB', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOB', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOB,datrc_4d)      
  else
    !print*,'RINOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOB,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOB', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPBO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPBO,datrc_4d)      
  else
    !print*,'RIPBO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPBO,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPBO', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROQCS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,ROQCS,datrc_4d)      
  else
    !print*,'ROQCS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROQCS,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROQCS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROQAS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,ROQAS,datrc_4d)      
  else
    !print*,'ROQAS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROQAS,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROQAS', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RINHOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHOR,datrc_3d)      
  else
    !print*,'RINHOR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHOR,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RINHOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RINOOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOOR,datrc_3d)      
  else
    !print*,'RINOOR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOOR,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RINOOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RIPOOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPOOR,datrc_3d)      
  else
    !print*,'RIPOOR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPOOR,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrbHetetrophCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RIPOOR', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMEhetr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element',units='g d-2', interpinic_flag='skip', &
      data=datpr5, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OMEhetr,datrc_5d)      
  else
    !print*,'OMC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMEhetr,datrc_5d)        
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMEhetr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  
  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROXYSff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,ROXYSff,datrc_3d)      
  else
    !print*,'ROXYSff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROXYSff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='ROXYSff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX4ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX4ff,datrc_3d)      
  else
    !print*,'RVMX4ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX4ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX4ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX3ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX3ff,datrc_3d)      
  else
    !print*,'RVMX3ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX3ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX3ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX2ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX2ff,datrc_3d)      
  else
    !print*,'RVMX2ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX2ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX2ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX1ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMX1ff,datrc_3d)      
  else
    !print*,'RVMX1ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMX1ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMX1ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB4ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB4ff,datrc_3d)      
  else
    !print*,'RVMB4ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB4ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB4ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB3ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB3ff,datrc_3d)      
  else
    !print*,'RVMB3ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB3ff,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB3ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB2ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMB2ff,datrc_3d)     
  else
    !print*,'RVMB2ff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMB2ff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMB2ff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHOff,datrc_3d)     
  else
    !print*,'RINHOff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHOff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOOff,datrc_3d)     
  else
    !print*,'RINOOff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOOff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPOOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPOOff,datrc_3d)     
  else
    !print*,'RIPOOff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPOOff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPOOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHBff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHBff,datrc_3d)     
  else
    !print*,'RINHBff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHBff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINHBff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOBff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOBff,datrc_3d)     
  else
    !print*,'RINOBff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOBff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RINOBff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPBOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPBOff,datrc_3d)     
  else
    !print*,'RIPBOff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPBOff,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutotrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RIPBOff', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)    
    call restartvar(ncid, flag, varname='RINHORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RINHORff,datrc_2d)     
  else
    !print*,'RINHORff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINHORff,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)       
    call restartvar(ncid, flag, varname='RINHORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)       
    call restartvar(ncid, flag, varname='RINOORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RINOORff,datrc_2d)     
  else
    !print*,'RINOORff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RINOORff,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)       
    call restartvar(ncid, flag, varname='RINOORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)       
    call restartvar(ncid, flag, varname='RIPOORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RIPOORff,datrc_2d)     
  else
    !print*,'RIPOORff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIPOORff,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutotrophCmplx)       
    call restartvar(ncid, flag, varname='RIPOORff', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)       
    call restartvar(ncid, flag, varname='OMEauto', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,OMEauto,datrc_4d)     
  else
    !print*,'OMCff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMEauto,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)           
    call restartvar(ncid, flag, varname='OMEauto', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='ORM', dim1name='column',dim2name='elements',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,ORM,datrc_5d)     
  else
    !print*,'ORC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ORM,datrc_5d)       
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='ORM', dim1name='column',dim2name='elements',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)               
    call restartvar(ncid, flag, varname='DOM', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic C micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM,datrc_4d)     
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM', dim1name='column',dim2name='ndoms',&
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP,datrc_4d)
  else
    !print*,'OQCH'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MacP', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter macropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='OHM', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OHM,datrc_4d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OHM,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='OHM', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      

  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='OSM', dim1name='column',dim2name='elements',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OSM,datrc_5d)
  else
    !print*,'OSC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OSM,datrc_5d)  
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='OSM', dim1name='column',dim2name='elements',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='OSA', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OSA,datrc_4d)
  else
    !print*,'OSA'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OSA,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='OSA', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  


  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMXC', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RVMXC,datrc_2d)
  else
    !print*,'RVMXC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RVMXC,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RVMXC', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ORGC', dim1name='column',dim2name='levsoi1',&
      long_name='total soil organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,ORGC,datrc_2d)
  else
    !print*,'ORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ORGC,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ORGC', dim1name='column',dim2name='levsoi1',&
      long_name='total soil organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ORGR', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,ORGR,datrc_2d)
  else
    !print*,'ORGR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ORGR,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ORGR', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_soil', dim1name='column',dim2name='fertN',&
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_soil,datrc_3d)
  else
    !print*,'FertN_soil'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_soil,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_soil', dim1name='column',dim2name='fertN', &
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_band', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_band,datrc_3d)
  else
    !print*,'FertN_band'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_band,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_band', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,WDNHB,datrc_2d)
  else
    !print*,'WDNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WDNHB,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNHB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,DPNHB,datrc_2d)
  else
    !print*,'DPNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPNHB,datrc_2d)  
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
    call cpcol(flag,NHW,NHE,NVN,NVS,WDNOB,datrc_2d)
  else
    !print*,'WDNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WDNOB,datrc_2d)  
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DPNOB,datrc_2d)
  else
    !print*,'DPNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPNOB,datrc_2d)  
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
    call cpcol(flag,NHW,NHE,NVN,NVS,WDPOB,datrc_2d)
  else
    !print*,'WDPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WDPOB,datrc_2d)  
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
    call cpcol(flag,NHW,NHE,NVN,NVS,DPPOB,datrc_2d)
  else
     !print*,'DPPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DPPOB,datrc_2d)  
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

  END SUBROUTINE wouts


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
    use ecosim_log_mod   , only : errMsg => shr_log_errMsg
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
    end if

    ! Branch run: 
    ! Restart file pathname is obtained from namelist "restartFileFullPath"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)

    if (nsrest==nsrBranch) then
       length = len_trim(restartFileFullPath)
       if (restartFileFullPath(length-2:length) == '.nc') then
          path = trim(restartFileFullPath) 
       else
          path = trim(restartFileFullPath) // '.nc'
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

  use EcoSIMCtrlMod     , only : etimer
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
       write(nio,'(a)')curr_date
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

  ! Do error checking on file
  
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
  call check_dim(ncid,'cansecz',NumOfLeafZenithSectors)
  call check_dim(ncid,'fertN',trc_confs%NFertNitro)
  call check_dim(ncid,'fertNb',trc_confs%NFertNitrob)
  call check_dim(ncid,'automicb',NumMicrobAutotrophCmplx)
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
  call check_dim(ncid,'hetrmicb',NumMicrbHetetrophCmplx)
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
  call ncd_defdim(ncid,'automicb',NumMicrobAutotrophCmplx, dimid)
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
  call ncd_defdim(ncid, 'sdim',3,dimid)   !grid dimension
  call ncd_defdim(ncid,'hetrmicb',NumMicrbHetetrophCmplx,dimid)
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
  nio = getavu()
  filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
  call opnfil( filename, nio, 'f' )  
  read(nio,*) fnamer
  read(nio,*)curr_date
  call relavu( nio )

  end subroutine get_restart_date
end module restartMod
