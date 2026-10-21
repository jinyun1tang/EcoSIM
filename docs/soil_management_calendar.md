# Soil-management calendar

The `&ecosim` option `soil_mgmt_use_model_year` defaults to `.false.`. This
preserves legacy tillage/fertilizer/irrigation schedule selection by the climate
data year. No existing namelist needs to change.

Set `soil_mgmt_use_model_year=.true.` to select the schedule NetCDF `year`
coordinate using the simulation year instead. Climate can still recycle its
forcing years without repeating a transient management event during spin-up.
The schedule must contain every requested model year, with explicit `NO` routes
for unmanaged years. Missing coordinates remain fatal. Supply event dates in
their actual model year, including the correct leap calendar.

This changes schedule selection only. Dose units, application calculations,
annual array resets and physical equations are unchanged. Fire, PFT management,
weather and GHG calendars are unaffected. It does not certify automatic
irrigation triggers, mixed automatic/scheduled operation or restart equivalence.
The selected option is recorded in the setup report.

`tests/test_management_reader_probe.py` compiles the actual management reader
and dispatcher against NetCDF, stubbing application consumers to check routing.
`tests/test_management_calendar_run.py` accepts an existing canonical CA-Let
native validation directory and a compiled executable. It runs three cold days
across 2001/2002 using recycled 2001 weather and a synthetic 0.01 g N/m2 event
only in the 2002 management coordinate. It checks recorded amendment mass and
exact default-versus-explicit-false history equality. This is a diagnostic test,
not a proposed site-management experiment or scientific calibration.

The artificial midyear cold test explicitly handles the driver's year-origin
stop counter and first history-bound initialization; the first cumulative
record is not interpreted as a daily flux. GESP's production cold candidates
require January-1 starts to avoid that separate calendar ambiguity. No model
clock or restart code is changed by the scheduling option.
