"""
Plot WCT dv/v (time x frequency) for all station pairs found in an MSNoise
output folder, several panels per page, as many pages as needed.

Pairs are discovered automatically from the .nc files, so you don't have to
list stations. Set `stations` to restrict to a subset.
"""
import matplotlib
matplotlib.use("Agg")  # no display needed (fixes the X11 "port 6012" error)

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import DateFormatter, YearLocator, MonthLocator
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path


class Config:
    def __init__(self):
        self.volcano = "Kilauea"
        self.component = "ZZ"
        self.window = "10d_1d"
        self.filter_id = "02"
        self.base_path = Path("WCT_MERGED")

        # Empty list = use every pair found in the folder.
        # Otherwise only pairs where BOTH stations are in the list.
        self.stations = []
        self.include_autocorr = False  # include STA-STA files if present

        # Time period
        self.startdate = "2022-01-01"
        self.enddate = "2026-12-31"
        self.titleref = "Static ref"

        # Layout: panels per page = nrows * ncols (4x2 = 8, 3x2 = 6)
        self.nrows = 4
        self.ncols = 2

        # Colour scale
        self.use_fixed_dvv_scale = False
        self.dvv_scale = {"min": -0.2, "max": 0.2}
        self.percentile = 98          # symmetric scale from |dv/v| percentile
        self.same_scale_all_pages = True

        # Coherence masking
        self.apply_coherence_mask = False
        self.coherence_threshold = 0.0

        # Events (vertical dashed lines)
        self.event_list = []          # e.g. [pd.to_datetime("2023-01-05")]

        # Output
        self.out_dir = Path("Figures")
        self.save_png = True
        self.save_pdf = True          # one multi-page PDF with all pages


def discover_pairs(cfg):
    folder = cfg.base_path / cfg.filter_id / cfg.window / cfg.component
    if not folder.exists():
        raise FileNotFoundError(f"Folder not found: {folder.resolve()}")

    pairs = []
    for f in sorted(folder.glob("*.nc")):
        try:
            id1, id2 = f.stem.split("_")
            sta1, sta2 = id1.split(".")[1], id2.split(".")[1]
        except (ValueError, IndexError):
            print(f"Skipping unrecognised file name: {f.name}")
            continue
        if sta1 == sta2 and not cfg.include_autocorr:
            continue
        if cfg.stations and not (sta1 in cfg.stations and sta2 in cfg.stations):
            continue
        pairs.append((f"{sta1}-{sta2}", f))

    if not pairs:
        raise RuntimeError(f"No matching .nc files in {folder.resolve()}")
    print(f"Found {len(pairs)} pairs in {folder}")
    return pairs


def load_pair(path, cfg):
    with xr.open_dataset(path) as ds:
        dvv = ds["dvv"].to_pandas()
        coh = ds["coh"].to_pandas() if "coh" in ds else None

    # Make sure rows are time, columns are frequency
    if not isinstance(dvv.index, pd.DatetimeIndex) and isinstance(dvv.columns, pd.DatetimeIndex):
        dvv = dvv.T
        coh = coh.T if coh is not None else None

    dvv.index = pd.to_datetime(dvv.index)
    dvv = dvv.sort_index()

    if cfg.apply_coherence_mask and coh is not None:
        coh.index = pd.to_datetime(coh.index)
        coh = coh.reindex_like(dvv)
        dvv = dvv.where(coh > cfg.coherence_threshold)

    dvv = dvv.loc[cfg.startdate:cfg.enddate]
    dvv.columns = dvv.columns.astype(float)
    dvv = dvv.sort_index(axis=1)

    if dvv.empty or dvv.isna().all().all():
        return None

    # One row per day so data gaps show up as gaps
    dvv.index = dvv.index.normalize()
    dvv = dvv.groupby(level=0).mean()
    full_index = pd.date_range(cfg.startdate, cfg.enddate, freq="D")
    return dvv.reindex(full_index)


def symmetric_limit(frames, pct):
    vals = [d.values.ravel() for d in frames if d is not None]
    if not vals:
        return 0.2
    m = np.nanpercentile(np.abs(np.concatenate(vals)), pct)
    return m if np.isfinite(m) and m > 0 else 0.2


def plot_page(items, page_no, n_pages, vlim, cfg):
    nr, nc = cfg.nrows, cfg.ncols
    fig, axs = plt.subplots(nr, nc, figsize=(10 * nc, 4.2 * nr),
                            sharex=True, squeeze=False)
    axs = axs.ravel()

    start, end = pd.to_datetime(cfg.startdate), pd.to_datetime(cfg.enddate)
    cmap = mpl.cm.RdBu_r
    norm = mpl.colors.Normalize(vmin=-vlim, vmax=vlim)

    for i, (ax, (label, dvv)) in enumerate(zip(axs, items)):
        ax.set_title(label, fontsize=16)
        if dvv is None:
            ax.text(0.5, 0.5, "No data in period", ha="center", va="center",
                    transform=ax.transAxes, fontsize=13)
        else:
            ax.pcolormesh(dvv.index, dvv.columns.values,
                          np.ma.masked_invalid(dvv.values.T),
                          cmap=cmap, norm=norm, shading="nearest")

        for ev in cfg.event_list:
            if start <= ev <= end:
                ax.axvline(ev, color="k", ls="--", lw=2)

        ax.set_xlim(start, end)
        ax.xaxis.set_major_locator(YearLocator())
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_major_formatter(DateFormatter("%Y"))
        ax.tick_params(axis="x", which="major", length=8, width=1.5,
                       labelsize=12, labelbottom=True)
        ax.tick_params(axis="x", which="minor", length=4)
        ax.tick_params(axis="y", labelsize=12)
        if i % nc == 0:
            ax.set_ylabel("Frequency (Hz)", fontsize=14)

    # Hide empty panels on the last page
    for ax in axs[len(items):]:
        ax.set_visible(False)

    fig.subplots_adjust(left=0.06, right=0.9, top=0.93, bottom=0.05,
                        hspace=0.35, wspace=0.12)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cbar.set_label("dv/v (%)", fontsize=16)
    cbar.ax.tick_params(labelsize=12)

    fig.suptitle(f"{cfg.volcano} | {cfg.component} | {cfg.window} | filter {cfg.filter_id}"
                 f" | {cfg.titleref} | page {page_no}/{n_pages}", fontsize=18)
    return fig


def main():
    cfg = Config()
    pairs = discover_pairs(cfg)

    data = []
    for label, path in pairs:
        print(f"Loading {label}")
        data.append((label, load_pair(path, cfg)))

    per_page = cfg.nrows * cfg.ncols
    pages = [data[i:i + per_page] for i in range(0, len(data), per_page)]

    if cfg.use_fixed_dvv_scale:
        global_lim = max(abs(cfg.dvv_scale["min"]), abs(cfg.dvv_scale["max"]))
    else:
        global_lim = symmetric_limit([d for _, d in data], cfg.percentile)

    cfg.out_dir.mkdir(exist_ok=True)
    stem = f"{cfg.volcano}_dvv_{cfg.component}_{cfg.window}_f{cfg.filter_id}"
    pdf = PdfPages(cfg.out_dir / f"{stem}.pdf") if cfg.save_pdf else None

    for p, items in enumerate(pages, start=1):
        if cfg.use_fixed_dvv_scale or cfg.same_scale_all_pages:
            vlim = global_lim
        else:
            vlim = symmetric_limit([d for _, d in items], cfg.percentile)

        fig = plot_page(items, p, len(pages), vlim, cfg)
        if cfg.save_png:
            fig.savefig(cfg.out_dir / f"{stem}_p{p:02d}.png", dpi=200)
        if pdf:
            pdf.savefig(fig)
        plt.close(fig)
        print(f"Page {p}/{len(pages)} done")

    if pdf:
        pdf.close()
    print(f"Output in {cfg.out_dir.resolve()}")


if __name__ == "__main__":
    main()
