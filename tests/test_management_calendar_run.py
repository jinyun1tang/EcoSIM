"""Cold three-day year-boundary integration test using supplied native inputs.

Temporary synthetic N application, not a scientific site-management scenario.
Climate year 2001 repeats while simulation advances from 2001 into 2002.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from netCDF4 import Dataset


def control(text, name, value):
    result, count = re.subn(rf'(?im)^\s*{name}\s*=.*$', f' {name}={value}', text)
    if count != 1:
        raise ValueError(f'one {name} required')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, template = args.output.resolve(), args.template.resolve()
    if root.exists():
        raise FileExistsError('fresh output required')
    root.mkdir(parents=True)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    report = {'binary_sha256': digest(args.executable), 'cases': {}}
    base = (template / 'CA-Let_validation.namelist').read_text()
    inputs = [name for name in re.findall(r"(?im)^\s*\w+_in\s*=\s*'([^']+)'", base) if name != 'NO']
    snapshots = {}
    for role, flag in (('default', None), ('false', '.false.'), ('model_year', '.true.')):
        folder = root / role
        folder.mkdir()
        for name in inputs:
            (folder / name).symlink_to(template / name)
        path = folder / 'calendar.nc'
        with Dataset(path, 'w') as ds:
            for name, size in (('year', 2), ('ntopou', 1), ('nchar1', 10), ('nfert', 12), ('ncharf', 128)):
                ds.createDimension(name, size)
            ds.createVariable('year', 'i4', ('year',))[:] = [2001, 2002]
            for name in ('NH1', 'NH2', 'NV1', 'NV2'):
                ds.createVariable(name, 'i4', ('ntopou',))[:] = 1
            for name in ('fertf', 'tillf', 'irrigf'):
                var = ds.createVariable(name, 'S1', ('year', 'ntopou', 'nchar1'))
                var[:] = np.array([[list('NO'.ljust(10))]] * 2, dtype='S1')
            ds['fertf'][1, 0] = np.array(list('fert2002'.ljust(10)), dtype='S1')
            tape = ds.createVariable('fert2002', 'S1', ('nfert', 'ncharf'))
            tape[:] = np.full((12, 128), b' ', dtype='S1')
            # NH4 dose = 0.01 g N/m2, all other doses/depth/codes zero.
            record = '01012002 0.01 ' + ' '.join(['0'] * 25)
            tape[0] = np.array(list(record.ljust(128)), dtype='S1')
        text = base
        for name, value in (('case_name', repr('calendar_' + role)), ('start_date', "'20011231000000'"),
                ('forc_periods', '2001,2001,3'), ('soil_mgmt_in', "'calendar.nc'"),
                # EcoSIM's cold clock offset counts from January 1. Day 367
                # is three simulated days after the December-31 start here.
                ('plant_model', '.false.'), ('stop_option', "'ndays'"), ('stop_n', '367'),
                ('rest_opt', "'never'"), ('hist_fincl1', "'AMENDED_N_col:I','ECO_GPP_col','tNH4_col'")):
            text = control(text, name, value)
        if flag is not None:
            text = text.replace("soil_mgmt_in='calendar.nc'", "soil_mgmt_in='calendar.nc'\n soil_mgmt_use_model_year=" + flag)
        (folder / 'case.namelist').write_text(text)
        checksum = digest(path)
        with (folder / 'ecosim.log').open('w') as log:
            result = subprocess.run([str(args.executable.resolve()), 'case.namelist'], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        log = (folder / 'ecosim.log').read_text(errors='replace')
        if result.returncode or re.search(r'ENDRUN:|Program received signal|error test failure', log, re.I):
            raise ValueError(f'{role}: native run failed; inspect its log')
        steps = set(map(int, re.findall(r'\bnstep\s*=\s*(\d+)', log, re.I)))
        # This driver logs steps at checkpoints, not at every hourly update.
        # Verify the absolute final step AND complete daily time bounds below.
        assert max(steps) == 8808, (min(steps), max(steps), len(steps))
        assert digest(path) == checksum
        arrays = {}
        for history in sorted(folder.glob('*.ecosim.h*.nc')):
            with Dataset(history) as ds:
                for name, var in ds.variables.items():
                    if np.issubdtype(var.dtype, np.number):
                        data = np.ma.asarray(var[:])
                        assert np.isfinite(data.compressed()).all(), name
                        arrays.setdefault(name, []).append(np.asarray(data))
        arrays = {key: np.concatenate(value) for key, value in arrays.items()}
        # HistFileMod initializes the first bound to zero, even for this
        # artificial late-year cold start. Do not interpret that first
        # cumulative-state record as a 365-day daily flux. The two subsequent
        # Jan-1/Jan-2 intervals must each cover exactly one day.
        assert np.array_equal(arrays['time_bounds'], [[0, 365], [365, 366], [366, 367]])
        amended = arrays['AMENDED_N_col'].reshape(-1)
        assert len(amended) == 3
        expected = [0, .01, .01] if flag == '.true.' else [0, 0, 0]
        assert np.allclose(amended, expected, rtol=1.e-5, atol=1.e-9), (role, amended.tolist())
        snapshots[role] = arrays
        report['cases'][role] = {'status': 'pass', 'hourly_steps': 72, 'amended_N_g_m2': amended.tolist(),
                                'history_time_bounds': arrays['time_bounds'].tolist()}
    assert snapshots['default'].keys() == snapshots['false'].keys()
    assert all(np.array_equal(snapshots['default'][name], snapshots['false'][name]) for name in snapshots['default'])
    report.update(status='pass', legacy_default_matches_false_exactly=True,
                  cold_start_note='First late-year cold record has an initialization bound from zero; not a daily flux. GESP production cold candidates require January 1.',
                  scope='cold year-boundary routing and small prescribed N dose; no restart or irrigation certification')
    (root / 'calendar_integration_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
