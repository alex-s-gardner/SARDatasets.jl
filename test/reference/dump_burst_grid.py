"""Dump the burst placements hyp3-autorift's merge produces, as golden values.

`burst_grid` reproduces the azimuth arithmetic of `merge_bursts_in_swath` and
`get_azimuth_reference_offsets` in hyp3-autorift's `s1_isce3.py`, so what it computes is checked
against what those functions compute. The two are transcribed here rather than imported: running them
needs ISCE3, COMPASS and GDAL, and the arithmetic under test is integer index arithmetic over
annotation fields, which needs none of that.

Run against a directory holding the granules named in `sentinel1_metadata*.json`:

    SLCDATASETS_S1_DIR=<dir> python3 test/reference/dump_burst_grid.py

Writes `burst_grid.json` beside this file.

Two departures from the reference are deliberate and recorded here rather than smoothed over.

The reference writes each burst into the merged array in turn, and where an overlap spans an odd
number of rows its floor division leaves one row that two consecutive bursts both write. The later
write survives. So `merge_end` is the extent a burst is *written* into, and `eff_merge_end` is the
extent it still *owns* afterwards; the effective extents tile the image exactly, and they are what a
reader placing each burst once must reproduce.

The reference then takes range samples with `slice(first_valid_sample, last_valid_sample)`, whose end
is exclusive, while `last_valid_sample` is an inclusive index — so it drops the last valid sample of
every burst. The values here are the inclusive bounds; the Julia reader keeps that sample, and its
tests assert the one-sample difference rather than inheriting it.
"""

import json
import os
import sys
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

HERE = os.path.dirname(os.path.abspath(__file__))
GOLDEN = ["sentinel1_metadata.json", "sentinel1_metadata_s1b.json"]


def annotation(zf, swath, polarization):
    """The product annotation of one subswath, from the zipped SAFE."""
    names = [
        n
        for n in zf.namelist()
        if "/annotation/" in n
        and n.endswith(".xml")
        and "calibration" not in n
        and "rfi" not in n
        and f"-iw{swath}-slc-{polarization}-" in n.split("/")[-1]
    ]
    if len(names) != 1:
        raise SystemExit(f"expected one IW{swath} {polarization} annotation, found {names}")
    return ET.fromstring(zf.read(names[0]))


def read_bursts(root):
    """Each burst's start and valid region, by s1reader's `burst_from_xml` reduction."""
    lines_per_burst = int(root.findtext("swathTiming/linesPerBurst"))
    samples_per_burst = int(root.findtext("swathTiming/samplesPerBurst"))
    interval = float(root.findtext("imageAnnotation/imageInformation/azimuthTimeInterval"))

    bursts = []
    for element in root.findall("swathTiming/burstList/burst"):
        first_samples = [int(v) for v in element.findtext("firstValidSample").split()]
        last_samples = [int(v) for v in element.findtext("lastValidSample").split()]
        first_valid_line = [v >= 0 for v in first_samples].index(True)
        n_valid_lines = [v >= 0 for v in first_samples].count(True)
        last_valid_line = first_valid_line + n_valid_lines - 1
        bursts.append(
            {
                "sensing_start": datetime.fromisoformat(element.findtext("azimuthTime")),
                "first_valid_line": first_valid_line,
                "last_valid_line": last_valid_line,
                "first_valid_sample": max(
                    first_samples[first_valid_line], first_samples[last_valid_line]
                ),
                "last_valid_sample": min(
                    last_samples[first_valid_line], last_samples[last_valid_line]
                ),
            }
        )
    return bursts, lines_per_burst, samples_per_burst, interval


def azimuth_reference_offsets(bursts, interval):
    """`get_azimuth_reference_offsets`: each burst's valid region as rows of the merged grid."""
    offsets = []
    sensing_start = bursts[0]["sensing_start"]
    for burst in bursts:
        start = burst["sensing_start"] + timedelta(seconds=burst["first_valid_line"] * interval)
        first = round((start - sensing_start).total_seconds() / interval)
        last = first + (burst["last_valid_line"] - burst["first_valid_line"]) + 1
        offsets.append([int(first), int(last)])
    return offsets


def merge(bursts, lines_per_burst, interval):
    """`merge_bursts_in_swath`: the merged height, and the slices each burst is written with."""
    n = len(bursts)
    burst_length = timedelta(seconds=(lines_per_burst - 1.0) * interval)
    sensing_start = bursts[0]["sensing_start"]
    sensing_end = bursts[-1]["sensing_start"] + burst_length
    num_az_lines = 1 + round((sensing_end - sensing_start).total_seconds() / interval)

    limits = azimuth_reference_offsets(bursts, interval)
    rows = []
    for i, burst in enumerate(bursts):
        limit = limits[i]
        head = 0 if i == 0 else (limits[i - 1][1] - limit[0]) // 2
        tail = 0 if i == n - 1 else (limit[1] - limits[i + 1][0]) // 2
        rows.append(
            {
                "burst": i + 1,
                "merge_start": int(limit[0] + head),
                "merge_end": int(limit[1] - tail),
                "burst_start": int(burst["first_valid_line"] + head),
                "burst_end": int(1 + burst["last_valid_line"] - tail),
                "first_valid_sample": burst["first_valid_sample"],
                "last_valid_sample": burst["last_valid_sample"],
            }
        )

    # What each burst still owns once the later write of a shared row has happened.
    for i, row in enumerate(rows):
        end = row["merge_end"] if i == n - 1 else min(row["merge_end"], rows[i + 1]["merge_start"])
        row["eff_merge_end"] = end
        row["eff_burst_end"] = row["burst_end"] - (row["merge_end"] - end)

    # The effective extents must tile the merged image without gap or overlap; anything else means
    # this transcription has drifted from the reference.
    for i in range(n - 1):
        if rows[i]["eff_merge_end"] != rows[i + 1]["merge_start"]:
            raise SystemExit(f"burst {i + 1} and {i + 2} do not meet")
    for row in rows:
        if row["eff_merge_end"] - row["merge_start"] != row["eff_burst_end"] - row["burst_start"]:
            raise SystemExit(f"burst {row['burst']} takes a different count than it gives")

    return {
        "num_az_lines": int(num_az_lines),
        "lines_per_burst": lines_per_burst,
        "rows": rows,
    }


def main():
    directory = os.environ.get("SLCDATASETS_S1_DIR", "")
    if not directory:
        raise SystemExit("set SLCDATASETS_S1_DIR to the directory holding the granules")

    out = {}
    for name in GOLDEN:
        with open(os.path.join(HERE, name)) as f:
            gold = json.load(f)
        safe = os.path.join(directory, gold["safe"])
        polarization = gold["polarization"]
        with zipfile.ZipFile(safe) as zf:
            per_swath = {}
            for swath in (1, 2, 3):
                root = annotation(zf, swath, polarization)
                bursts, lines, samples, interval = read_bursts(root)
                merged = merge(bursts, lines, interval)
                merged["samples_per_burst"] = samples
                per_swath[str(swath)] = merged
        out[gold["safe"]] = per_swath

    dest = os.path.join(HERE, "burst_grid.json")
    with open(dest, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"wrote {dest}")
    for safe, swaths in out.items():
        heights = {s: v["num_az_lines"] for s, v in swaths.items()}
        stacked = {s: v["lines_per_burst"] * len(v["rows"]) for s, v in swaths.items()}
        print(f"  {safe}")
        print(f"    merged   {heights}")
        print(f"    stacked  {stacked}")


if __name__ == "__main__":
    sys.exit(main())
