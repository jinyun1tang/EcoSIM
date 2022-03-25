  SUBROUTINE splitc(NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS)
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils, only : endrun
  use GridConsts
  use EcoSIMHistMod
  use GridDataType

  implicit none
  integer, intent(in) :: NT,NE,NAX,NDX,NTX,NEX,NHW,NHE,NVN,NVS

  character(len=1024) :: str
  integer :: nz,nx,ny,n
  integer :: failure
  character(len=*), parameter :: modfile='splitc'
! execution begins here

  nz=1
  do nx=nhw,nhe
  do ny=nvn,nvs
    if(nz .lt. np0(ny,nx))nz=np0(ny,nx)
  enddo
  enddo
  do N=1,10
	  if(datac(N+20,NE,NEX) .NE. 'NO')then
	    close((N+30))
      close((N+40))
!     call splits(NHW,NHE,NVN,NVS,OUTS(N))
      call splits(NHW,NHE,NVN,NVS,outdir,OUTS(N), failure)
      if(failure==1)call endrun('Fail to process file '// &
        trim(outdir)//trim(OUTS(N))//' in '//trim(modfile))
      str='rm -f ' // OUTS(N)
      call system (str)
!     call splitp(NHW,NHE,NVN,NVS,nz,OUTP(N))
      call splitp(NHW,NHE,NVN,NVS,nz,outdir, OUTP(N), failure)
      if(failure==1)call endrun('Fail to process file '// &
        trim(outdir)//trim(OUTP(N))//' in '//trim(modfile))
      str = 'rm -f ' // OUTP(N)
      call system (str)
    endif
  enddo
  RETURN
  END subroutine splitc
