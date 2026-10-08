"""
Broadband albedo for several stars and surfaces, in two bands and in six bands.

Based on Broadband_albedo_calculator_2_band.py and
broadband_albedo_calculator_6_band.py (Vidya Venkatesan), which are python
versions of ice_gcm.pro. The albedo calculation is the same; this script runs
it for every star/surface pair, prints the band albedos and saves figures:
  two_band_albedo.png  stellar spectra in wavelength bins above the surface
                       reflectance and its two-band albedo for each star
  six_band_albedo.png  the same for the six bands
  weighted_albedo.png  the flux-weighted interpolated albedo for each pair

For each pair the script
  1. finds the overlapping wavelength range of the two files,
  2. interpolates both spectra onto a common 10000-point grid
     (linear for the stellar flux, cubic spline for the reflectance),
  3. integrates reflectance x flux over each band and divides by the
     integrated flux in that band.
The two bands are split at 0.7 microns and cover the whole overlapping range.
The six bands are the surface-albedo bands of ROCKE-3D, VIS and NIR1 to NIR5,
from 0.3 to 4.0 microns. A six-band value that the files only
partly cover is marked with *, and one they do not cover at all is given as
n/a (the original six-band script prints 0 there).

Usage:
    python plot_weighted_albedo.py                 (two bands and six bands)
    python plot_weighted_albedo.py --bands 2       (two bands only)
    python plot_weighted_albedo.py --bands 6       (six bands only)
    python plot_weighted_albedo.py --data-dir /path/to/spectra
    python plot_weighted_albedo.py --no-show       (save the figures, open no window)

Inputs (two columns each, extra columns and header lines are ignored):
    stellar file = wavelength, flux (any units, the flux is normalized)
    surface file = wavelength, reflectance (0 - 1)
Wavelength units are detected from the smallest wavelength in each file and
converted to microns. Set "unit" in the tables below to 'm', 'um', 'nm' or 'A'
if the detection reported on screen is wrong for one of your files.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import scipy.interpolate as interpol
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter

# --------------------------------------------------------------------------
# What to run. Edit these two tables to change the stars or the surfaces.
# --------------------------------------------------------------------------
STARS = [
    # file,                  label,               unit,  color,     linestyle
    ("hd128167_scaled.txt", "F2V HD128167",       None, "#2f9fd0", "-"),
    ("sun_scaled.txt",      "G2V Sun",            None, "#dba400", "-"),
    ("hd22049_scaled.txt",  "K2V HD22049",        None, "#e8601c", "-"),
]

SURFACES = [
    # file,                       label,                                unit,  color,     linestyle
    ("snow_bluemarine_50_50.txt", "50% snow / 50% blue marine ice",     None, "#1a1a19", "-"),
    ("CO2_i200.txt",              "CO$_2$ ice, 200 $\\mu$m grains",     None, "#5b6477", "--"),
]

BANDS = "both"     # "2", "6" or "both": which band sets to print and plot
BAND_SPLIT = 0.7   # microns, boundary between the two bands

# The six surface-albedo bands of ROCKE-3D (Way et al. 2017, Table 3): one
# visible band and five near-infrared bands. That table starts the VIS band at
# 0.33 microns and the text of the paper at 0.30; this keeps the 0.30 of the
# six-band calculator. Change the first edge to 0.33 to follow the table.
SIX_BAND_EDGES = [0.3, 0.77, 0.86, 1.25, 1.5, 2.2, 4.0]   # microns
SIX_BAND_NAMES = ["VIS", "NIR1", "NIR2", "NIR3", "NIR4", "NIR5"]
NGRID = 10000      # points in the common wavelength grid
SNORM = 1360       # W m^-2, total flux each stellar spectrum is scaled to

# Two-band and six-band figures
X_RANGE = (0.2, 5.0)    # microns, wavelength range of the figure (log axis)
BINS_PER_DECADE = 20    # wavelength bins the stellar spectra are averaged into

# Weighted-albedo figure
X_MAX = None       # microns, right edge of the surface and weighted panels
                   # (None = full overlap range)
STAR_X_MAX = 5.0   # microns, right edge of the stellar spectra panel

FORTRAN_EXPONENT = re.compile(r"(\d\.?)[dD]([+-]?\d)")
TO_MICRONS = {"m": 1.0e6, "um": 1.0, "nm": 1.0e-3, "A": 1.0e-4}

# The band sets. The two-band set has open ends: it uses the whole range that
# the star and surface files share, as the two-band calculator does.
SCHEMES = {
    "2": {
        "name": "two-band",
        "model": "",
        "edges": [-np.inf, BAND_SPLIT, np.inf],
        "headers": ["VIS", "NIR"],
        "ranges": [f"< {BAND_SPLIT:g} $\\mu$m", f"$\\geq$ {BAND_SPLIT:g} $\\mu$m"],
        "plain": [f"< {BAND_SPLIT:g} um", f">= {BAND_SPLIT:g} um"],
    },
    "6": {
        "name": "six-band",
        "model": " (ROCKE-3D bands)",
        "edges": SIX_BAND_EDGES,
        "headers": SIX_BAND_NAMES,
        "ranges": [f"{lo:.2f}\u2013{hi:.2f}" for lo, hi in zip(SIX_BAND_EDGES[:-1], SIX_BAND_EDGES[1:])],
        "plain": [f"{lo:.2f}-{hi:.2f}" for lo, hi in zip(SIX_BAND_EDGES[:-1], SIX_BAND_EDGES[1:])],
    },
}


def read_text(path):
    """Read a text file whatever its encoding (UTF-8, UTF-16 or Latin-1)."""
    raw = Path(path).read_bytes()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in raw[:200]:
        return raw.decode("utf-16", errors="replace")
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def to_float(text):
    """float() that also accepts Fortran exponents such as 1.5D-07."""
    return float(FORTRAN_EXPONENT.sub(r"\1E\2", text))


def read_two_columns(path):
    """Read the first two numeric columns of a text file, skipping headers."""
    x, y = [], []
    lines = read_text(path).splitlines()
    for line in lines:
        parts = line.replace(",", " ").replace(";", " ").split()
        if len(parts) < 2:
            continue
        try:
            a, b = to_float(parts[0]), to_float(parts[1])
        except ValueError:
            continue  # header or comment line
        if np.isfinite(a) and np.isfinite(b):
            x.append(a)
            y.append(b)
    if len(x) < 2:
        shown = [line[:100] + ("..." if len(line) > 100 else "") for line in lines[:5]]
        preview = "\n".join(f"    {text!r}" for text in shown) or "    (the file is empty)"
        raise ValueError(
            f"{path}: found fewer than two rows of the form 'wavelength flux'.\n"
            f"  The file is {Path(path).stat().st_size} bytes and has {len(lines)} lines. "
            f"It starts with:\n{preview}")
    x, y = np.array(x), np.array(y)
    # The interpolation needs strictly increasing wavelengths.
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    keep = np.concatenate(([True], np.diff(x) > 0))
    return x[keep], y[keep]


def guess_unit(wavelength):
    """Guess the wavelength unit from the smallest wavelength in the file."""
    shortest = np.min(wavelength[wavelength > 0])
    if shortest < 1.0e-3:
        return "m"
    if shortest < 10.0:
        return "um"
    if shortest < 1000.0:
        return "nm"
    return "A"


def load_spectrum(data_dir, filename, unit):
    """Return wavelength in microns and the second column of the file."""
    path = Path(data_dir) / filename
    wavelength, values = read_two_columns(path)
    detected = unit is None
    if detected:
        unit = guess_unit(wavelength)
    wavelength = wavelength * TO_MICRONS[unit]
    how = "detected" if detected else "set"
    print(f"  {filename}: {len(wavelength)} points, wavelength unit {unit} ({how}), "
          f"{wavelength.min():.3f} - {wavelength.max():.3f} microns")
    return wavelength, values


def broadband_albedo(lamda, flux, wave, albedo):
    """The 2-band calculation for one star and one surface."""
    # Overlapping wavelength range of the stellar and albedo data
    start_g = max(np.min(wave), np.min(lamda))
    end_g = min(np.max(wave), np.max(lamda))
    if end_g <= start_g:
        raise ValueError("the stellar and surface spectra do not overlap in wavelength")

    # Common wavelength grid
    n = np.arange(NGRID)
    wavelengthgrid = start_g + (end_g - start_g) * n / (NGRID - 1)
    dlamda = (np.max(wavelengthgrid) - np.min(wavelengthgrid)) / (len(wavelengthgrid) - 1)

    # Interpolate both spectra onto the grid
    stellarInterpolate = np.interp(wavelengthgrid, lamda, flux)
    albedoInterpolate = interpol.CubicSpline(wave, albedo)(wavelengthgrid)

    # Scale the stellar flux so that it integrates to SNORM
    total_sed = np.sum(stellarInterpolate * dlamda)
    new_flux = stellarInterpolate * (SNORM / total_sed)

    def band(mask):
        total = np.sum(new_flux[mask] * dlamda)
        if total == 0:
            return 0.0
        return np.sum(albedoInterpolate[mask] * new_flux[mask] * dlamda) / total

    return {
        "grid": wavelengthgrid,
        "albedo_interp": albedoInterpolate,
        # Fraction of the stellar flux per micron: integrates to 1 over the grid
        "weight": new_flux / SNORM,
        "start": start_g,
        "end": end_g,
        "band": band,
        "lower": band(wavelengthgrid < BAND_SPLIT),
        "upper": band(wavelengthgrid >= BAND_SPLIT),
        "total": band(np.ones(NGRID, dtype=bool)),
    }


def style_axis(ax):
    ax.grid(True, color="#e4e3df", linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#8a8984")
    ax.tick_params(colors="#52514e", labelsize=9)


def mark_band_split(ax, label=False):
    ax.axvline(BAND_SPLIT, color="#8a8984", linewidth=1.0, linestyle=":")
    if label:
        ax.annotate(f"{BAND_SPLIT} $\\mu$m band split", xy=(BAND_SPLIT, 1.0),
                    xycoords=("data", "axes fraction"), xytext=(4, -4),
                    textcoords="offset points", va="top", ha="left",
                    fontsize=8, color="#52514e")


def band_albedos(r, edges):
    """Albedo in each band between consecutive edges, and over all of them.

    r is the result of broadband_albedo for one star and one surface.
    """
    grid = r["grid"]

    def one(lo, hi):
        mask = (grid >= lo) & (grid < hi)
        span_lo, span_hi = max(lo, r["start"]), min(hi, r["end"])
        covered = span_hi > span_lo and mask.any()
        cut = ((np.isfinite(lo) and span_lo > lo + 1e-9) or
               (np.isfinite(hi) and span_hi < hi - 1e-9))
        return {"lo": span_lo, "hi": span_hi,
                "value": r["band"](mask) if covered else np.nan,
                "partial": bool(covered and cut)}

    return {"bands": [one(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])],
            "total": one(edges[0], edges[-1])}


def cell(band, digits=3):
    """Text for one albedo value: n/a if not covered, * if partly covered."""
    if np.isnan(band["value"]):
        return "n/a"
    return f"{band['value']:.{digits}f}" + ("*" if band["partial"] else "")


def print_bands(stars, surfaces, results, key):
    scheme = SCHEMES[key]
    plain = lambda text: text.replace("$", "").replace("\\mu", "u").replace("_", "")
    width = max(12, max(len(name) for name in scheme["plain"]) + 3)
    print(f"{scheme['name'].capitalize()} albedo{scheme['model']}")
    flagged = False
    for a_label, *_ in surfaces:
        print(f"  {plain(a_label)}")
        print(f"    {'star':<18}" + "".join(f"{name:>{width}}" for name in scheme["headers"])
              + f"{'total':>{width}}   files overlap (um)")
        print(f"    {'':<18}" + "".join(f"{name:>{width}}" for name in scheme["plain"]))
        for s_label, *_ in stars:
            r = results[(s_label, a_label)]
            table = r[key]
            cells = [cell(b, 6) for b in table["bands"]] + [cell(table["total"], 6)]
            flagged = flagged or any("*" in c or c == "n/a" for c in cells)
            print(f"    {s_label:<18}" + "".join(f"{c:>{width}}" for c in cells)
                  + f"   {r['start']:.2f} - {r['end']:.2f}")
    if flagged:
        print("  * band only partly inside the range shared by the star and surface files;"
              " n/a: outside it")
    print()


def bin_edges(x_lo, x_hi, breaks):
    """Log-spaced bin edges that include every band edge in breaks."""
    points = [x_lo] + sorted(b for b in breaks if x_lo < b < x_hi) + [x_hi]
    edges = [np.array([x_lo])]
    for lo, hi in zip(points[:-1], points[1:]):
        n = max(1, int(round(BINS_PER_DECADE * np.log10(hi / lo))))
        edges.append(np.logspace(np.log10(lo), np.log10(hi), n + 1)[1:])
    return np.concatenate(edges)


def scaled_cumulative_flux(lamda, flux):
    """Running integral of the flux, scaled to SNORM over the whole file."""
    running = np.concatenate(([0.0], np.cumsum(0.5 * (flux[1:] + flux[:-1]) * np.diff(lamda))))
    return running * SNORM / running[-1]


def band_figure(stars, surfaces, results, key):
    """Binned stellar spectra above the surface reflectance and band albedos."""
    scheme = SCHEMES[key]
    x_lo, x_hi = X_RANGE
    shade = np.clip(scheme["edges"], x_lo, x_hi)
    breaks = [e for e in scheme["edges"] if x_lo < e < x_hi]
    two_band = key == "2"

    # Table layout: surfaces side by side if the columns fit, otherwise stacked
    n_cols = len(scheme["headers"]) + 1
    stacked = n_cols * len(surfaces) > 8
    block_rows = len(stars) + 3   # title, band names, wavelength ranges, stars
    flagged = any(b["partial"] or np.isnan(b["value"])
                  for r in results.values() for b in r[key]["bands"] + [r[key]["total"]])
    table_rows = (len(surfaces) * block_rows + 0.6 * (len(surfaces) - 1)) if stacked else block_rows
    table_rows += 1.2 if flagged else 0
    table_height = 0.25 * table_rows

    fig_height = 9.3 + table_height
    fig = plt.figure(figsize=(8.2, fig_height))
    # Rows: spectra, albedo, a gap for the wavelength label, the table
    grid = fig.add_gridspec(4, 1, height_ratios=[3.6, 3.6, 0.35, table_height], hspace=0.12)
    ax_star = fig.add_subplot(grid[0])
    ax_alb = fig.add_subplot(grid[1], sharex=ax_star)
    ax_tab = fig.add_subplot(grid[3])

    tints = ("#f5f3ec", "#e4e9f3")
    for ax in (ax_star, ax_alb):
        for i, (lo, hi) in enumerate(zip(shade[:-1], shade[1:])):
            if hi > lo:
                ax.axvspan(lo, hi, color=tints[i % 2], zorder=0)
        for edge in breaks:
            ax.axvline(edge, color="#52514e", linewidth=0.9 if two_band else 0.6,
                       linestyle="-.", zorder=1)
        ax.grid(True, axis="y" if not two_band else "both", color="#cfcec8", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(colors="#52514e", labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("#8a8984")
    for name, lo, hi in zip(scheme["headers"], shade[:-1], shade[1:]):
        if hi > lo:
            ax_star.text(np.sqrt(lo * hi), 1.03, name, transform=ax_star.get_xaxis_transform(),
                         ha="center", va="bottom", fontsize=12 if two_band else 10,
                         fontweight="bold", color="#1a1a19")
    # Top panel: each spectrum scaled to SNORM over its whole file, then averaged
    # in wavelength bins so that every bin keeps its share of the flux.
    edges = bin_edges(x_lo, x_hi, breaks)
    for label, color, ls, lamda, flux in stars:
        running = scaled_cumulative_flux(lamda, flux)
        binned = np.diff(np.interp(edges, lamda, running)) / np.diff(edges)
        covered = (edges[:-1] >= lamda.min()) & (edges[1:] <= lamda.max())
        binned = np.where(covered, binned, np.nan)
        if two_band:
            label = f"{label}   {np.interp(BAND_SPLIT, lamda, running) / SNORM:.0%}"
        ax_star.step(edges, np.append(binned, binned[-1]), where="post", color=color,
                     linewidth=1.8, label=label)
    ax_star.set_ylabel("Stellar flux at the planet (W m$^{-2}$ $\\mu$m$^{-1}$)", fontsize=10, fontweight="bold")
    ax_star.set_ylim(bottom=0)
    ax_star.legend(frameon=not two_band, framealpha=0.85, edgecolor="none", fontsize=9,
                   loc="upper right", alignment="left", title_fontsize=9,
                   title=f"Share of starlight below {BAND_SPLIT:g} $\\mu$m" if two_band else None)
    plt.setp(ax_star.get_xticklabels(), visible=False)

    # Bottom panel: reflectance spectrum of each surface, and on top of it the
    # band albedos that each star gives for that surface.
    first_star = stars[0][0]
    handles = []
    for a_label, a_color, a_ls, _, _ in surfaces:
        r = results[(first_star, a_label)]
        ax_alb.plot(r["grid"], r["albedo_interp"], color=a_color, linestyle=a_ls, linewidth=1.1)
        handles.append(Line2D([], [], color=a_color, linestyle=a_ls, linewidth=1.1,
                              label=f"{a_label}: reflectance spectrum"))
        for s_label, s_color, _, _, _ in stars:
            x, y = [], []
            for band in results[(s_label, a_label)][key]["bands"]:
                missing = np.isnan(band["value"])
                x += [band["lo"], band["hi"]] if not missing else [np.nan, np.nan]
                y += [band["value"], band["value"]]
            ax_alb.plot(x, y, color=s_color, linestyle=a_ls, linewidth=2.2)
    handles.append(Line2D([], [], color="#8a8984", linewidth=2.2,
                          label=f"Thick steps: {scheme['name']} albedo for each star\n"
                                "(star's colour, surface's line style)"))
    ax_alb.legend(handles=handles, frameon=True, framealpha=0.85, edgecolor="none",
                  fontsize=8.5, loc="best")
    ax_alb.set_ylabel("Surface albedo", fontsize=10, fontweight="bold")
    ax_alb.set_xlabel("Wavelength ($\\mu$m)", fontsize=10, fontweight="bold")
    ax_alb.set_ylim(0, 1)

    ax_alb.set_xscale("log")
    ax_alb.set_xlim(x_lo, x_hi)
    if two_band:
        ticks = {0.2, 0.3, 0.5, 1, 2, 3, 5, 10, BAND_SPLIT}
    else:
        ticks = {x_lo, x_hi, *scheme["edges"]}
    ax_alb.xaxis.set_major_locator(FixedLocator(sorted(t for t in ticks if x_lo <= t <= x_hi)))
    ax_alb.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    ax_alb.xaxis.set_minor_formatter(NullFormatter())
    if not two_band:
        ax_alb.tick_params(axis="x", which="minor", bottom=False)
        plt.setp(ax_alb.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=8.5)

    # Table of the flux-weighted albedos. Rows count down from the top.
    ax_tab.axis("off")
    ax_tab.set_xlim(0, 1)
    ax_tab.set_ylim(table_rows, 0)
    ink, muted = "#1a1a19", "#52514e"
    names = scheme["headers"] + ["Total"]
    ranges = scheme["ranges"] + [""]
    if not stacked:
        ax_tab.text(0.0, 0.5, "Flux-weighted albedo", fontsize=9, fontweight="bold", va="center", color=ink)
    for j, (a_label, *_) in enumerate(surfaces):
        if stacked:
            top, left, width = j * (block_rows + 0.6), 0.26, 0.74
            ax_tab.text(0.0, top + 0.5, f"Flux-weighted albedo: {a_label}", fontsize=9,
                        fontweight="bold", va="center", color=ink)
            ax_tab.text(0.035, top + 2.5, "Wavelength ($\\mu$m)", fontsize=8.5, va="center", color=muted)
        else:
            top, width = 0, 0.74 / len(surfaces)
            left = 0.26 + j * width
            ax_tab.text(left + 0.5 * width, 0.5, a_label, fontsize=9, fontweight="bold",
                        ha="center", va="center", color=ink)
        for i, (s_label, s_color, *_) in enumerate(stars):
            if stacked or j == 0:
                ax_tab.plot([0.012], [top + i + 3.5], marker="s", markersize=7, color=s_color, clip_on=False)
                ax_tab.text(0.035, top + i + 3.5, s_label, fontsize=9, va="center", color=ink)
            table = results[(s_label, a_label)][key]
            for k, band in enumerate(table["bands"] + [table["total"]]):
                ax_tab.text(left + (k + 0.5) * width / len(names), top + i + 3.5, cell(band),
                            fontsize=9, ha="center", va="center", color=ink)
        for k, (name, wavelengths) in enumerate(zip(names, ranges)):
            x = left + (k + 0.5) * width / len(names)
            ax_tab.text(x, top + 1.5, name, fontsize=9, fontweight="bold", ha="center", va="center", color=ink)
            ax_tab.text(x, top + 2.5, wavelengths, fontsize=8.5, ha="center", va="center", color=muted)
    if flagged:
        ax_tab.text(0.0, table_rows - 0.4, "* band only partly inside the wavelength range shared by the "
                    "star and surface files;  n/a: outside it", fontsize=8, va="center", color=muted)

    fig.suptitle(f"Stellar spectra and the {scheme['name']} surface albedo{scheme['model']}", fontsize=13,
                 fontweight="bold", y=1 - 0.3 / fig_height)
    fig.subplots_adjust(left=0.11, right=0.97, top=1 - 0.95 / fig_height, bottom=0.2 / fig_height)
    return fig


def weighted_figure(stars, surfaces, results):
    """Stellar spectra, surface reflectance and the flux-weighted albedo."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    (ax_star, ax_surf), weighted_axes = axes
    x_min = min(r["start"] for r in results.values())
    x_max = max(r["end"] for r in results.values())

    # Stellar spectra in W m^-2 um^-1. Each one is scaled so that it integrates
    # to SNORM over the whole wavelength range of its file, which puts all the
    # stars on the same total flux and lets their peaks differ.
    for label, color, ls, lamda, flux in stars:
        total = np.sum(0.5 * (flux[1:] + flux[:-1]) * np.diff(lamda))
        shown = lamda <= STAR_X_MAX
        ax_star.plot(lamda[shown], flux[shown] * SNORM / total, color=color, linestyle=ls,
                     linewidth=1.3, label=label)
    ax_star.set_title(f"Stellar spectra, each scaled to {SNORM} W m$^{{-2}}$", fontsize=11, loc="left")
    ax_star.set_xlabel("Wavelength [$\\mu$m]", fontsize=10)
    ax_star.set_ylabel("Stellar flux [W/(m$^2$ $\\mu$m)]", fontsize=10)
    ax_star.set_ylim(bottom=0)
    ax_star.legend(frameon=False, fontsize=9)

    # Interpolated reflectance of each surface (the same for every star)
    first_star = stars[0][0]
    for label, color, ls, _, _ in surfaces:
        r = results[(first_star, label)]
        ax_surf.plot(r["grid"], r["albedo_interp"], color=color, linestyle=ls, linewidth=1.6, label=label)
    ax_surf.set_title("Interpolated surface reflectance", fontsize=11, loc="left")
    ax_surf.set_xlabel("Wavelength [$\\mu$m]", fontsize=10)
    ax_surf.set_ylabel("Reflectance", fontsize=10)
    ax_surf.set_ylim(bottom=0)
    ax_surf.legend(frameon=False, fontsize=9)

    # Flux-weighted albedo: one panel per surface, one line per star
    for ax, (a_label, *_) in zip(weighted_axes, surfaces):
        for s_label, color, ls, _, _ in stars:
            r = results[(s_label, a_label)]
            ax.plot(r["grid"], r["albedo_interp"] * r["weight"], color=color, linestyle=ls,
                    linewidth=1.6, label=f"{s_label}: albedo = {r['total']:.3f}")
        ax.set_title(f"Weighted albedo: {a_label}", fontsize=11, loc="left")
        ax.set_xlabel("Wavelength [$\\mu$m]", fontsize=10)
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False, fontsize=9, title="Area under curve = broadband albedo",
                  title_fontsize=8.5, alignment="left")
    weighted_axes[0].set_ylabel("Reflectance $\\times$ flux fraction ($\\mu$m$^{-1}$)", fontsize=10)
    top = max(ax.get_ylim()[1] for ax in weighted_axes)
    for ax in weighted_axes:
        ax.set_ylim(0, top)

    for i, ax in enumerate(axes.flat):
        style_axis(ax)
        mark_band_split(ax, label=(i == 0))
        ax.set_xlim(x_min, x_max if X_MAX is None else X_MAX)
    ax_star.set_xlim(0, STAR_X_MAX)

    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0].strip())
    parser.add_argument("--data-dir", default=Path(__file__).resolve().parent,
                        help="folder with the spectra (default: the folder of this script)")
    parser.add_argument("--bands", choices=["2", "6", "both"], default=BANDS,
                        help=f"band sets to print and plot (default: {BANDS})")
    parser.add_argument("--out", default="two_band_albedo.png",
                        help="file for the two-band figure (default: two_band_albedo.png)")
    parser.add_argument("--out-six", default="six_band_albedo.png",
                        help="file for the six-band figure (default: six_band_albedo.png)")
    parser.add_argument("--out-weighted", default="weighted_albedo.png",
                        help="file for the weighted-albedo figure (default: weighted_albedo.png)")
    parser.add_argument("--no-show", action="store_true",
                        help="save the figure without opening a window")
    args = parser.parse_args()

    print("Reading stellar spectra")
    stars = [(label, color, ls) + load_spectrum(args.data_dir, name, unit)
             for name, label, unit, color, ls in STARS]
    print("Reading surface reflectance spectra")
    surfaces = []
    for name, label, unit, color, ls in SURFACES:
        wave, albedo = load_spectrum(args.data_dir, name, unit)
        if np.max(albedo) > 1.5:
            print(f"    reflectance goes up to {np.max(albedo):.1f}: read as percent, divided by 100")
            albedo = albedo / 100.0
        surfaces.append((label, color, ls, wave, albedo))

    # Run the calculation for every star/surface pair
    keys = ["2", "6"] if args.bands == "both" else [args.bands]
    results = {}
    for s_label, _, _, lamda, flux in stars:
        for a_label, _, _, wave, albedo in surfaces:
            r = broadband_albedo(lamda, flux, wave, albedo)
            for key in keys:
                r[key] = band_albedos(r, SCHEMES[key]["edges"])
            results[(s_label, a_label)] = r

    print()
    for key in keys:
        print_bands(stars, surfaces, results, key)

    figures = [(band_figure(stars, surfaces, results, key), {"2": args.out, "6": args.out_six}[key])
               for key in keys]
    figures.append((weighted_figure(stars, surfaces, results), args.out_weighted))
    for fig, name in figures:
        out = Path(name)
        fig.savefig(out, dpi=200)
        print(f"Figure saved to {out.resolve()}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()