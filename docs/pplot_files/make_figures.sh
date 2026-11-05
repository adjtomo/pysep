#!/usr/bin/env bash
# Regenerate the example figures in the PrettyPlot docs page (docs/pplot.md)
#
# Usage: run from a directory holding the example MiniSEED day files and
# StationXML (see "Example data" in docs/pplot.md), e.g.
#
#     cd /path/to/example_data
#     bash /path/to/pysep/docs/pplot_files/make_figures.sh
#
# Needs the `pplot` command from PySEP. If it is not on your PATH, set PPLOT,
# e.g. PPLOT=~/micromamba/envs/pysep/bin/pplot
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${HERE}"
PPLOT="${PPLOT:-pplot}"
COMMON=(--noshow --dpi 150)

run() {
    # Run PrettyPlot with the common options, quietly unless it fails
    echo "-> $1"
    local name="$1"; shift
    MPLBACKEND=Agg "${PPLOT}" "$@" "${COMMON[@]}" \
        --save "${OUT}/${name}.png" > /dev/null
}

EVENT=(--xlim 2025-03-31T12:00:00 2025-03-31T12:03:00)
WINDY=(--xlim 2025-04-08T07:30:00 2025-04-08T10:30:00)

# 1. Raw data, relative time axis in seconds from the start of the file
run 01_raw XX.MITDC..HHZ.2025.090 --xlim 43200 43380

# 2. Response removal, bandpass, absolute time in local time (HST)
run 02_response XX.MITDC..HHZ.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 1 --fmax 20 --time a-10 \
    --xlim 2025-03-31T02:00:00 2025-03-31T02:03:00

# 3. Several stations overlain
run 03_overlay XX.MITD[CUW]..HHZ.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 1 --fmax 20 --time a \
    --xlim 2025-03-31T12:00:20 2025-03-31T12:00:50 --colors 'C?' \
    --labels Covered Uncovered Windscreen --linewidth 0.75

# 4. Record section of the three components, normalized
run 04_recsec XX.MITDC..HH?.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 1 --fmax 20 --time a \
    "${EVENT[@]}" --wf_type recsec --wf_norm --wf_spacing_exact \
    --wf_recsec_spacing 2 --colors C0 C1 C2 --ylim -1.2 6.5 \
    --ylabel "Normalized velocity"

# 5. Grid of subplots, one per file, each with its own y-limits
run 05_grid XX.MITDC..HH?.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 1 --fmax 20 --time a \
    "${EVENT[@]}" --wf_type grid --colors C0 C1 C2 --fig_size 7 7

# 6. Spectrogram of a 3-hour window as the wind drops
run 06_spectrogram XX.MITDC..HHZ.2025.098 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 0.5 --time a "${WINDY[@]}" \
    --spectrogram --sp_vmax -130 --sp_clip 0.5 1

# 7. Log-frequency spectrogram of the event with a custom window length
run 07_spectrogram_log XX.MITDC..HHZ.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --time a "${EVENT[@]}" \
    --spectrogram --sp_logscale --sp_wlen 8 --sp_cmap magma \
    --sp_clip 0.4 1

# 8. Wind speed from the co-located anemometer on a second y-axis
run 08_add_trace XX.MITDC..HHZ.2025.098 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 0.5 --time a "${WINDY[@]}" \
    --add_trace XX.MITDA..LWS.2025.098 --tr_ylabel "Wind speed (m/s)" \
    --tr_label "MITDA wind speed" --tr_color C3 --ylim 2E-5

# 9. Time marks
run 09_tmarks XX.MITDC..HHZ.2025.090 --inv ISHTAR.xml \
    --rr_pre_filt 0.05 0.1 40 45 --fmin 1 --fmax 20 --time a \
    "${EVENT[@]}" --tmarks 2025-03-31T12:00:28 2025-03-31T12:00:45 \
    --tmarks_c r b

# 10. Merge a gappy file (two short gaps near 02:10:25 UTC)
run 10_merge HV.MITD..HHZ.2025.090 --merge --time a \
    --xlim 2025-03-31T02:09:30 2025-03-31T02:11:30 \
    --tmarks 2025-03-31T02:10:22.96 2025-03-31T02:10:35.33 --tmarks_c r

echo "figures written to ${OUT}"
