"""Compile a real NetCDF conversion probe against the exact I/O guard source."""
import argparse
import json
import shlex
import signal
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default="gfortran")
    parser.add_argument("--nf-config", default="nf-config")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    def flags(option):
        return shlex.split(subprocess.check_output([args.nf_config, option], text=True))
    with tempfile.TemporaryDirectory(prefix="ecosim-output-precision-") as temporary:
        work = Path(temporary)
        executable = work / "probe"
        command = [args.fc, "-O0", "-g", "-std=f2008", "-ffree-line-length-none",
                   "-ffpe-trap=invalid,zero,overflow,underflow", *flags("--fflags"),
                   str(root / "f90src/IOutils/HistoryOutputPrecisionMod.F90"),
                   str(root / "tests/history_output_precision_probe.F90"), *flags("--flibs"), "-o", str(executable)]
        subprocess.run(command, cwd=work, check=True)
        results = {}
        for mode in ("guarded", "existing-flag", "unguarded", "physics", "overflow"):
            folder = work / mode
            folder.mkdir()
            outcome = subprocess.run([str(executable), mode], cwd=folder, capture_output=True, text=True, check=False)
            passed = (outcome.returncode == 0 and "PASS:" in outcome.stdout if mode in ("guarded", "existing-flag")
                      else outcome.returncode in (-signal.SIGFPE, -signal.SIGILL))
            results[mode] = {"pass": passed, "returncode": outcome.returncode,
                             "stdout": outcome.stdout, "stderr": outcome.stderr}
        report = {"status": "pass" if all(row["pass"] for row in results.values()) else "fail", "checks": results}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
