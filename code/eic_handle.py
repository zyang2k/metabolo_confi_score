import numpy as np
import matplotlib.pyplot as plt
from pyteomics import mzml


def extract_eic(mzml_path, target_mz, tolerance=0.005, tolerance_unit="da",
                ms_level=1, rt_range=None):
    """Extract the EIC trace for a target m/z from an mzML file.

    Returns
    -------
    rts : np.ndarray
        Retention times in minutes.
    intensities : np.ndarray
        Summed intensity within the m/z window for each scan.
    """
    if tolerance_unit.lower() == "ppm":
        delta = target_mz * tolerance / 1e6
    elif tolerance_unit.lower() == "da":
        delta = tolerance
    else:
        raise ValueError("tolerance_unit must be 'ppm' or 'da'")

    mz_low, mz_high = target_mz - delta, target_mz + delta

    rts, intensities = [], []
    with mzml.read(mzml_path) as reader:
        for spectrum in reader:
            if spectrum.get("ms level") != ms_level:
                continue

            scan_info = spectrum["scanList"]["scan"][0]
            rt = scan_info["scan start time"]
            if scan_info["scan start time"].unit_info == "second":
                rt = float(rt) / 60.0
            else:
                rt = float(rt)

            if rt_range is not None and not (rt_range[0] <= rt <= rt_range[1]):
                continue

            mz_array = spectrum["m/z array"]
            int_array = spectrum["intensity array"]
            mask = (mz_array >= mz_low) & (mz_array <= mz_high)
            rts.append(rt)
            intensities.append(int_array[mask].sum() if mask.any() else 0.0)

    return np.asarray(rts), np.asarray(intensities)


def plot_eic(mzml_path, target_mz, tolerance=0.005, tolerance_unit="da",
             ms_level=1, rt_range=None, ax=None, show=True, interactive=False):
    """Read an mzML file and plot the extracted ion chromatogram for a given m/z.

    Parameters
    ----------
    mzml_path : str
        Path to the mzML file.
    target_mz : float
        Target m/z value to extract.
    tolerance : float
        m/z extraction window half-width.
    tolerance_unit : {"ppm", "da"}
        Units for the tolerance value.
    ms_level : int
        MS level to extract from (default MS1).
    rt_range : tuple of float, optional
        (rt_min, rt_max) in minutes to restrict the chromatogram.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on (matplotlib path only). Ignored when interactive=True.
    show : bool
        Render the figure (plt.show / fig.show).
    interactive : bool
        If True, render an interactive plotly figure with hover tooltips.

    Returns
    -------
    If interactive=False (default):
        (rts, intensities) tuple of np.ndarray.
    If interactive=True:
        (rts, intensities, plotly.graph_objects.Figure).
    """
    rts, intensities = extract_eic(
        mzml_path, target_mz, tolerance=tolerance,
        tolerance_unit=tolerance_unit, ms_level=ms_level, rt_range=rt_range,
    )

    title = f"EIC m/z {target_mz:.4f} ± {tolerance} {tolerance_unit}"
    has_signal = len(intensities) and intensities.max() > 0
    if has_signal:
        apex_idx = int(np.argmax(intensities))
        apex_rt = float(rts[apex_idx])
        apex_int = float(intensities[apex_idx])

    if interactive:
        import plotly.graph_objects as go

        hover = [
            f"RT: {r:.3f} min ({r*60:.1f} s)<br>Intensity: {i:.3e}"
            for r, i in zip(rts, intensities)
        ]
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=rts, y=intensities, mode="lines",
            line=dict(width=1.5),
            hovertext=hover, hoverinfo="text",
            name="EIC",
        ))
        if has_signal:
            fig.add_trace(go.Scatter(
                x=[apex_rt], y=[apex_int], mode="markers+text",
                marker=dict(symbol="triangle-down", color="red", size=12),
                text=[f"{apex_rt*60:.1f} s<br>{apex_int:.2e}"],
                textposition="bottom right",
                hovertext=[f"APEX<br>RT: {apex_rt:.3f} min ({apex_rt*60:.1f} s)"
                           f"<br>Intensity: {apex_int:.3e}"],
                hoverinfo="text",
                name="apex",
            ))
        fig.update_layout(
            title=title,
            xaxis_title="Retention time (min)",
            yaxis_title="Intensity",
            hovermode="x unified",
            template="simple_white",
            width=900, height=450,
        )
        if show:
            fig.show()
        return rts, intensities, fig

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(rts, intensities, lw=1)
    ax.set_xlabel("Retention time (min)")
    ax.set_ylabel("Intensity")
    ax.set_title(title)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    if has_signal:
        ax.plot(apex_rt, apex_int, "rv", markersize=8, zorder=5)
        ax.annotate(f"{apex_rt * 60:.1f} s\n{apex_int:.2e}",
                    xy=(apex_rt, apex_int),
                    xytext=(5, -5), textcoords="offset points",
                    va="top", ha="left", fontsize=9, color="red")

    if show:
        plt.show()

    return rts, intensities
