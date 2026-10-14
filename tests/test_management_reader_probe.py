"""Compile actual reader excerpts against NetCDF; stub only event consumers.

The routing test executes ReadManagementFiles verbatim, so the real Fortran
NetCDF indexing is exercised. It is not an irrigation-process validation.
"""
import argparse
import json
import re
import shlex
import subprocess
import tempfile
from pathlib import Path

PRELUDE = """
module EcoSIMCtrlMod
  character(len=256) :: soil_mgmt_in='management.nc'
  logical :: Lirri_auto=.false.
  logical :: soil_mgmt_use_model_year=.false.
end module
module probe_reader
  use netcdf
  use EcoSIMCtrlMod
  implicit none
  integer,parameter :: ncd_nowrite=nf90_nowrite
  type :: file_desc_t
    integer :: fh
  end type
  type :: Var_desc_t
    integer :: varid
  end type
  logical :: first_topou=.false.,lverb=.false.
  character(len=*),parameter :: mod_filename='probe'
  character(len=10) :: seen(3,2)=''
contains
  subroutine endrun(message,line)
    character(len=*),intent(in)::message
    integer,intent(in),optional::line
    print *, trim(message)
    error stop 1
  end subroutine
  function int2str(value) result(output)
    integer,intent(in)::value
    character(len=32)::output
    write(output,'(i0)')value
  end function
  subroutine PrintInfo(message)
    character(len=*),intent(in)::message
  end subroutine
  subroutine check_ret(code,message)
    integer,intent(in)::code
    character(len=*),intent(in)::message
    if(code/=nf90_noerr)call endrun(nf90_strerror(code))
  end subroutine
  subroutine ncd_pio_openfile(file,path,mode)
    type(file_desc_t),intent(out)::file
    character(len=*),intent(in)::path
    integer,intent(in)::mode
    call check_ret(nf90_open(path,mode,file%fh),'open')
  end subroutine
  subroutine ncd_pio_closefile(file)
    type(file_desc_t),intent(inout)::file
    call check_ret(nf90_close(file%fh),'close')
  end subroutine
  function get_dim_len(file,name) result(length)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name
    integer::length,dim
    call check_ret(nf90_inq_dimid(file%fh,name,dim),'dimension')
    call check_ret(nf90_inquire_dimension(file%fh,dim,len=length),'length')
  end function
  subroutine ncd_getvar(file,name,index,value)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name
    integer,intent(in)::index
    integer,intent(out)::value
    integer::vid,data(1)
    call check_ret(nf90_inq_varid(file%fh,name,vid),'var')
    call check_ret(nf90_get_var(file%fh,vid,data,start=[index],count=[1]),'get')
    value=data(1)
  end subroutine
  subroutine check_var(file,name,description,found)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name
    type(Var_desc_t),intent(out)::description
    logical,intent(out)::found
    found=nf90_inq_varid(file%fh,name,description%varid)==nf90_noerr
  end subroutine
  subroutine ReadTillageFile(file,name,nh1,nh2,nv1,nv2,dimension)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name,dimension
    integer,intent(in)::nh1,nh2,nv1,nv2
    seen(2,nh1)=name
  end subroutine
  subroutine ReadFertlizerFile(file,name,nh1,nh2,nv1,nv2)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name
    integer,intent(in)::nh1,nh2,nv1,nv2
    seen(1,nh1)=name
  end subroutine
  subroutine ReadIrrigationFile(file,name,nh1,nh2,nv1,nv2)
    type(file_desc_t),intent(in)::file
    character(len=*),intent(in)::name
    integer,intent(in)::nh1,nh2,nv1,nv2
    seen(3,nh1)=name
  end subroutine
"""
EPILOGUE = """
end module
program management_probe
  use probe_reader
  implicit none
  character(len=16)::mode
  integer::unit
  call get_command_argument(1,mode)
  if(trim(mode)=='missing')then
    call ReadManagementFiles(2023)
    error stop 'missing year incorrectly accepted'
  endif
  first_topou=trim(mode)=='first'
  call ReadManagementFiles(2024)
  do unit=1,merge(1,2,first_topou)
    if(seen(1,unit)/='fert0'//trim(int2str(unit)))error stop 'fertilizer routing'
    if(seen(2,unit)/='till0'//trim(int2str(unit)))error stop 'tillage routing'
    if(seen(3,unit)/='irrig0'//trim(int2str(unit)))error stop 'irrigation routing'
  enddo
  if(first_topou.and.any(seen(:,2)/=''))error stop 'first_topou leaked'
  print *, 'PASS: native schedule routing'
end program
"""

CALENDAR_EPILOGUE = """
end module
program calendar_probe
  use probe_reader
  implicit none
  soil_mgmt_in='calendar.nc'
  ! Legacy mode deliberately reuses 2024 events in model year 2023.
  call ReadSoilManagementForYear(2023,2024)
  if(seen(1,1)/='fert01')error stop 'legacy default changed'
  seen=''
  soil_mgmt_use_model_year=.true.
  call ReadSoilManagementForYear(2023,2024)
  if(any(seen/=''))error stop 'transient event leaked into spinup'
  call ReadSoilManagementForYear(2024,2024)
  if(seen(1,1)/='fert01'.or.seen(1,2)/='fert02')error stop 'intended events missing'
  seen=''
  call ReadSoilManagementForYear(2025,2024)
  if(any(seen/=''))error stop 'transient event repeated with recycled climate'
  print *, 'PASS: model-year schedule isolated from recycled weather'
end program
"""

AUTO_PRELUDE = """
program auto_solutes_probe
  use iso_fortran_env,only: r8=>real64
  implicit none
  integer::I,IDY=0,NY=1,NX=1
  integer,parameter::idsalt_Al=1,idsalt_Fe=2,idsalt_Ca=3,idsalt_Mg=4,idsalt_Na=5
  integer,parameter::idsalt_K=6,idsalt_SO4=7,idsalt_Cl=8
  real(r8)::PHQ(366,1,1)=-1,NH4_irrig_mole_conc(366,1,1)=-1
  real(r8)::NO3_irrig_mole_conc(366,1,1)=-1,H2PO4_irrig_mole_conc(366,1,1)=-1
  real(r8)::trcsalt_irrig_mole_conc_col(8,366,1,1)=-1
  real(r8)::PHQX=6,CN4QX=28,CNOQX=42,CPOQX=62
  real(r8)::CALQX=54,CFEQX=111.6,CCAQX=80,CMGQX=48.6
  real(r8)::CNAQX=46,CKAQX=78.2,CSOQX=64,CCLQX=71
"""
AUTO_EPILOGUE = """
  if(any(PHQ/=6))error stop 'pH days not initialized'
  if(any(NH4_irrig_mole_conc/=2))error stop 'NH4 days not initialized'
  if(any(NO3_irrig_mole_conc/=3))error stop 'NO3 days not initialized'
  if(any(H2PO4_irrig_mole_conc/=2))error stop 'PO4 days not initialized'
  if(any(abs(trcsalt_irrig_mole_conc_col-2._r8)>1.e-6_r8))error stop 'salt days not initialized'
  print *, 'PASS: 366 daily solute entries initialized'
end program
"""


def main():
    import numpy as np
    from netCDF4 import Dataset

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--nf-config", default="nf-config")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = (root / "f90src/IOutils/ReadManagementMod.F90").read_text()
    routine = re.search(r"(?is)  subroutine ReadManagementFiles\(yeari\).*?end subroutine ReadManagementFiles", source)[0]
    selector = re.search(r"(?is)  subroutine ReadSoilManagementForYear\(.*?end subroutine ReadSoilManagementForYear", source)[0]
    with tempfile.TemporaryDirectory(prefix="ecosim-management-probe-") as directory:
        work = Path(directory)
        (work / "reader.F90").write_text(PRELUDE + routine + EPILOGUE)
        with Dataset(work / "management.nc", "w") as ds:
            for name, size in (("year", 1), ("ntopou", 2), ("nchar1", 10)):
                ds.createDimension(name, size)
            ds.createVariable("year", "i4", ("year",))[:] = [2024]
            for name in ("NH1", "NH2", "NV1", "NV2"):
                ds.createVariable(name, "i4", ("ntopou",))[:] = [1, 2] if name.startswith("NH") else [1, 1]
            for route, prefix in (("fertf", "fert"), ("tillf", "till"), ("irrigf", "irrig")):
                var = ds.createVariable(route, "S1", ("year", "ntopou", "nchar1"))
                var[:] = np.array([[list(f"{prefix}0{i}".ljust(10)) for i in (1, 2)]], dtype="S1")
        flags = lambda option: shlex.split(subprocess.check_output([args.nf_config, option], text=True))
        executable = work / "probe"
        subprocess.run([args.fc, "-O0", "-g", "-fcheck=all", "-ffree-line-length-none", *flags("--fflags"),
                        "reader.F90", *flags("--flibs"), "-o", str(executable)], cwd=work, check=True)
        checks = {}
        for mode in ("all", "first", "missing"):
            result = subprocess.run([str(executable), mode], cwd=work, capture_output=True, text=True, check=False)
            passed = (result.returncode != 0 and "management forcing year not present" in result.stdout
                      if mode == "missing" else result.returncode == 0 and "PASS:" in result.stdout)
            checks[mode] = {"pass": passed, "returncode": result.returncode}
        # Read the actual NetCDF year coordinate while exercising the actual
        # dispatcher. Stubbed consumers observe routing, not process response.
        with Dataset(work / "management.nc") as original, Dataset(work / "calendar.nc", "w") as ds:
            for name, dim in original.dimensions.items():
                ds.createDimension(name, 3 if name == "year" else len(dim))
            ds.createVariable("year", "i4", ("year",))[:] = [2023, 2024, 2025]
            for name, var in original.variables.items():
                if name == "year":
                    continue
                target = ds.createVariable(name, var.dtype, var.dimensions)
                if name in ("fertf", "irrigf", "tillf"):
                    target[:] = np.array([[list("NO".ljust(10))] * 2] * 3, dtype="S1")
                    target[1] = var[0]
                else:
                    target[:] = var[:]
        (work / "calendar.F90").write_text(PRELUDE + selector + "\n" + routine + CALENDAR_EPILOGUE)
        calendar_exe = work / "calendar_probe"
        subprocess.run([args.fc, "-O0", "-g", "-fcheck=all", "-ffree-line-length-none", *flags("--fflags"),
                        "calendar.F90", *flags("--flibs"), "-o", str(calendar_exe)], cwd=work, check=True)
        result = subprocess.run([str(calendar_exe)], cwd=work, capture_output=True, text=True, check=False)
        checks["model_year_with_recycled_weather"] = {"pass": result.returncode == 0 and "PASS:" in result.stdout,
                                                     "returncode": result.returncode}
        # Confirm the same spatial-index correction covers both fire tapes.
        fire = re.search(r"(?is)  subroutine ReadFire\(.*?end subroutine ReadFire", source)[0]
        checks["fire_indices"] = {"pass": len(re.findall(r"start\s*=\s*\(/1,ntopo/\)", fire)) == 2}
        loop = re.search(r"(?is)D220:\s*DO I = 1, 366.*?ENDDO D220", source)[0]
        (work / "auto.F90").write_text(AUTO_PRELUDE + loop + AUTO_EPILOGUE)
        auto = work / "auto_probe"
        subprocess.run([args.fc, "-O0", "-fcheck=all", "-ffree-line-length-none", "auto.F90", "-o", str(auto)],
                       cwd=work, check=True)
        result = subprocess.run([str(auto)], cwd=work, capture_output=True, text=True, check=False)
        checks["automatic_solute_initialization"] = {"pass": result.returncode == 0 and "PASS:" in result.stdout,
                                                      "returncode": result.returncode}
        report = {"status": "pass" if all(row["pass"] for row in checks.values()) else "fail", "checks": checks,
                  "scope": "actual schedule reader and NetCDF slicing; fire indexes source checked; not irrigation-process certification"}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
