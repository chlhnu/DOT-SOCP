from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional
import numpy as np
import matplotlib.pyplot as plt
import logging

from matplotlib.widgets import Slider
from matplotlib import animation
from dotsocp.solvers.utils.dot_barrier import _evaluate_barrier
from dotsocp.configs.icons import LOG_ICON
from matplotlib.colors import LinearSegmentedColormap, Normalize


_EVOLUTION_TYPES_1D = ["tile", "join"]
_EVOLUTION_TYPES_2D = ["imshow", "contourf", "contour", "contour3", "mesh"]


_parula_data = [
    [0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905], 
    [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143], 
    [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952, 
    0.779247619], [0.1252714286, 0.3242428571, 0.8302714286], 
    [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238, 
    0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571], 
    [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571, 
    0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429], 
    [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667, 
    0.8467], [0.0779428571, 0.5039857143, 0.8383714286], 
    [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571, 
    0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429], 
    [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524, 
    0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048, 
    0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667], 
    [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381, 
    0.7607190476], [0.0383714286, 0.6742714286, 0.743552381], 
    [0.0589714286, 0.6837571429, 0.7253857143], 
    [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429], 
    [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429, 
    0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048], 
    [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619, 
    0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667], 
    [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524, 
    0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905], 
    [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476, 
    0.4493904762], [0.609852381, 0.7473142857, 0.4336857143], 
    [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333], 
    [0.7184095238, 0.7411333333, 0.3904761905], 
    [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667, 
    0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762], 
    [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217], 
    [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857, 
    0.2886428571], [0.9738952381, 0.7313952381, 0.266647619], 
    [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857, 
    0.2164142857], [0.9955333333, 0.7860571429, 0.196652381], 
    [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857], 
    [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309], 
    [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333, 
    0.0948380952], [0.9661, 0.9514428571, 0.0755333333], 
    [0.9763, 0.9831, 0.0538]
]

from matplotlib.colors import LinearSegmentedColormap
PARULA_CMAP = LinearSegmentedColormap.from_list('parula', _parula_data)


def get_valid_evolution_types(dim: int) -> list[str]:
    """
    Get valid evolution visualization types based on problem dimension.

    Parameters
    ----------
    dim : int
        Problem dimension (1 or 2)

    Returns
    -------
    list[str]
        List of valid visualization types
    """
    if dim == 1:
        return _EVOLUTION_TYPES_1D
    elif dim == 2:
        return _EVOLUTION_TYPES_2D
    else:
        raise ValueError("Dimension must be 1 or 2.")


def show_evolution_1d(rho, show_func: Optional[str] = None, fig_name: str = "Density Evolution Over Time"):
    """
    Show evolution of 1D density rho over time.
    
    Parameters
    ----------
    rho : ndarray
        2D array of shape (Nx, Nt) containing density values
    show_func : Optional[str]
        Visualization type. One of `_EVOLUTION_TYPES_1D`. Defaults to "join".
    fig_name : str, optional
        Figure window name
    """
    logging.info(f"{LOG_ICON.get('stop', '[history]')}  Plotting density evolution...")

    if show_func is None:
        show_func = "join"
    elif show_func not in _EVOLUTION_TYPES_1D:
        raise ValueError("Invalid show_func: must be one of " + ", ".join(_EVOLUTION_TYPES_1D))
    
    Nx, Nt = rho.shape
    hx = 1.0 / Nx
    xx = np.linspace(hx / 2, 1 - hx / 2, Nx)
    
    y_min = 0.0
    y_max = 1.1 * np.max(rho)
    
    if show_func == "tile":
        fig = plt.figure(figsize=(12, 8))
        num_tile_horizon = 3
        num_tile_vertical = 3
        
        time_indices = np.round(np.linspace(0, Nt - 1, num_tile_vertical * num_tile_horizon)).astype(int)
        
        for idx, t in enumerate(time_indices, 1):
            ax = plt.subplot(num_tile_vertical, num_tile_horizon, idx)
            ax.plot(xx, rho[:, t], linewidth=1.5)
            ax.set_title(f"$t = {t}$")
            ax.set_ylim((y_min, y_max))
        
        plt.tight_layout()
    
    elif show_func == "join":
        fig = plt.figure(figsize=(9, 6))
        ax = plt.gca()
        
        num_curves = min(Nt, 10)
        tt = np.round(np.linspace(0, Nt - 1, num_curves)).astype(int)
        
        cmap = plt.get_cmap('turbo')
        colors = [cmap(i / (num_curves - 1)) for i in range(num_curves)]
        
        for ind, t in enumerate(tt):
            ax.plot(xx, rho[:, t], color=colors[ind], linewidth=1.5)
        
        ax.set_title(fig_name)
        ax.set_xlabel("Spatial Range")
        labels = [f"$t={100 * t / (Nt - 1):.2f}(\\%)$" for t in tt]
        ax.legend(labels, loc='best')
        ax.set_ylim((y_min, y_max))
    
    _adjust_fig(fig)
    
    return fig


def _prepare_evolution_2d_plot(
    rho: np.ndarray,
    show_func: Optional[str],
    fig_name: str,
    barrier_fn: Callable[[np.ndarray], np.ndarray] | None,
    cmap: str | None,
    show_colorbar: bool = True,
) -> tuple[plt.Figure, Callable[[int], object], int, str]:
    logging.info(f"{LOG_ICON.get('stop', '[history]')}  Plotting density evolution...")

    if show_func is None:
        show_func = "imshow"
    elif show_func not in _EVOLUTION_TYPES_2D:
        raise ValueError(f"Invalid show_func: must be one of " + ", ".join(_EVOLUTION_TYPES_2D))

    rho = np.asarray(rho, dtype=np.float64)
    if rho.ndim != 3:
        raise ValueError("rho must be a 3D array shaped (ny, nx, nt)")

    ny, nx, nt = rho.shape
    x = np.linspace(0.0, 1.0, nx)
    y = np.linspace(0.0, 1.0, ny)
    xx, yy = np.meshgrid(x, y)

    if barrier_fn is not None:
        barrier_grid = _evaluate_barrier(barrier_fn, xx.T, yy.T).T
        barrier_grid = np.asarray(barrier_grid, dtype=np.float64)
        barrier_mask = np.repeat(barrier_grid[:, :, np.newaxis], nt, axis=2)
    else:
        barrier_mask = None

    if cmap is None:
        if show_func == "imshow":
            cmap = "Greys_r" if barrier_fn is None else "turbo"
        elif show_func in {"contourf", "contour"}:
            cmap = "turbo"
        elif show_func == "mesh":
            cmap = PARULA_CMAP
        else:
            cmap = "viridis"

    if cmap == "parula":
        _cmap = PARULA_CMAP
    else:
        _cmap = plt.get_cmap(cmap)

    fig = plt.figure(figsize=(6.0, 6.0))
    ax = fig.add_subplot(111, projection="3d") if show_func in {"contour3", "mesh"} else fig.add_subplot(111)

    max_val = float(np.nanmax(rho)) if np.isfinite(rho).any() else 1.0

    def _apply_barrier_slice(data_slice: np.ndarray, mask_slice: np.ndarray | None, sentinel: float) -> np.ndarray:
        if mask_slice is None:
            return data_slice
        return np.where(mask_slice, sentinel, data_slice)

    artist = None
    colorbar = None
    contour_levels = None
    rho_imshow = rho
    rho_scaled = None
    contour_norm = None
    contour_cmap = _cmap

    def _get_dynamic_norm(data_slice, mask_slice):
        """Compute dynamic normalization for current frame, mimicking MATLAB's CLim behavior."""
        if mask_slice is not None:
            valid_data = data_slice[~np.asarray(mask_slice, dtype=bool)]
        else:
            valid_data = data_slice
        
        if valid_data.size > 0:
            vmin, vmax = float(np.min(valid_data)), float(np.max(valid_data))
        else:
            vmin, vmax = 0.0, 1.0
        
        return Normalize(vmin=vmin, vmax=vmax + 1e-9 if vmax <= vmin else vmax)

    if show_func == "imshow":
        if barrier_mask is not None:
            rho_imshow = np.where(barrier_mask, np.nan, rho)
        artist = ax.imshow(rho_imshow[:, :, 0], origin="lower", vmin=0.0, vmax=max_val, cmap=_cmap)
        if show_colorbar:
            colorbar = fig.colorbar(artist, ax=ax)
            colorbar.ax.yaxis.set_ticks([])
    elif show_func in {"contourf", "contour"}:
        num_levels = 128 if show_func == "contourf" else 30
        contour_levels = np.exp(np.linspace(0.0, np.log(255.0), num_levels))
        rho_scaled = rho * (255.0 / max_val) if max_val > 0 else rho.copy()
        contour_cmap = _cmap

        if barrier_mask is not None:
            contour_levels = np.concatenate(([-10.0], contour_levels))

        mask0 = None if barrier_mask is None else barrier_mask[:, :, 0]
        contour_norm = _get_dynamic_norm(rho_scaled[:, :, 0], mask0)
        data0 = _apply_barrier_slice(rho_scaled[:, :, 0], mask0, -np.inf)
        
        if show_func == "contourf":
            artist = ax.contourf(xx, yy, data0, levels=contour_levels, cmap=contour_cmap, norm=contour_norm, extend='max')
            if show_colorbar:
                colorbar = fig.colorbar(artist, ax=ax)
                colorbar.ax.yaxis.set_ticks([])
        else:
            artist = ax.contour(xx, yy, data0, levels=contour_levels, cmap=contour_cmap, norm=contour_norm, extend='max')
            if show_colorbar:
                colorbar = fig.colorbar(artist, ax=ax)
                colorbar.ax.yaxis.set_ticks([])
    elif show_func == "contour3":
        ax.view_init(elev=35, azim=-35)
        data0 = _apply_barrier_slice(rho[:, :, 0], None if barrier_mask is None else barrier_mask[:, :, 0], max_val)
        artist = ax.contour3D(xx, yy, data0, 30, cmap=_cmap)
        ax.set_zlim(0.0, max_val)
    elif show_func == "mesh":
        artist = ax.plot_surface(
            xx, yy, rho[:, :, 0], cmap=_cmap, rstride=1, cstride=1, linewidth=0,
            edgecolor="none",antialiased=False,shade=False,
        )
        ax.set_zlim(0.0, max_val)


    ax.set_title(fig_name)
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
    )

    if show_func in {"contour3", "mesh"}:
        ax.view_init(elev=30, azim=-120, roll=0)
        ax.grid(False)
        ax.set_axis_off()


    def _clear_contour_artist(art):
        if art is None:
            return
        collections = getattr(art, "collections", None)
        if collections is not None:
            for coll in collections:
                coll.remove()
            return
        try:
            art.remove()
        except Exception:
            pass

    def update_frame(t_index: int):
        nonlocal artist
        if nt == 0:
            return artist
        t_index = int(np.clip(t_index, 0, nt - 1))
        mask_slice = None if barrier_mask is None else barrier_mask[:, :, t_index]
        if show_func == "imshow":
            artist.set_data(rho_imshow[:, :, t_index])
        elif show_func == "contourf":
            _clear_contour_artist(artist)
            data = _apply_barrier_slice(rho_scaled[:, :, t_index], mask_slice, -np.inf)
            current_norm = _get_dynamic_norm(rho_scaled[:, :, t_index], mask_slice)
            artist = ax.contourf(xx, yy, data, levels=contour_levels, cmap=contour_cmap, norm=current_norm, extend='max')
        elif show_func == "contour":
            _clear_contour_artist(artist)
            data = _apply_barrier_slice(rho_scaled[:, :, t_index], mask_slice, -np.inf)
            current_norm = _get_dynamic_norm(rho_scaled[:, :, t_index], mask_slice)
            artist = ax.contour(xx, yy, data, levels=contour_levels, cmap=contour_cmap, norm=current_norm, extend='max')
        elif show_func == "contour3":
            ax.clear()
            ax.view_init(elev=35, azim=-35)
            ax.set_zlim(0.0, max_val)
            data = _apply_barrier_slice(rho[:, :, t_index], mask_slice, max_val)
            artist = ax.contour3D(xx, yy, data, 30, cmap=_cmap)
        elif show_func == "mesh":
            ax.clear()
            ax.set_zlim(0.0, max_val)
            artist = ax.plot_surface(
                xx, yy, rho[:, :, t_index], cmap=_cmap, rstride=1, cstride=1, linewidth=0,
                edgecolor="none",antialiased=False,shade=False,
            )

        if show_func in {"contour3", "mesh"}:
            ax.grid(False)
            ax.set_axis_off()

        _adjust_fig(fig)
        return artist

    return fig, update_frame, nt, show_func


def show_evolution_2d(
    rho: np.ndarray,
    show_func: Optional[str] = None,
    fig_name: str = "Density evolution",
    barrier_fn: Callable[[np.ndarray], np.ndarray] = None,
    cmap: str = None,
):
    """
    Visualize a time series of 2D densities using several display modes.

    Parameters
    ----------
    rho : ndarray
        Density array shaped (ny, nx, nt)
    show_func : Optional[str]
        Visualization type. One of `_EVOLUTION_TYPES_2D`. Defaults to "imshow".
    fig_name : str
        Figure window name
    barrier_fn : Optional[Callable[[np.ndarray], np.ndarray]]
        Function returning a barrier mask
    cmap : Optional[str]
        Colormap for visualization. See matplotlib colormaps.
    """
    fig, update_frame, nt, show_func = _prepare_evolution_2d_plot(
        rho=rho,
        show_func=show_func,
        fig_name=fig_name,
        barrier_fn=barrier_fn,
        cmap=cmap,
    )

    # Interactive slider (keep reference on fig to avoid GC)
    frame_slider = None
    if nt > 1:
        if show_func == "imshow":
            fig.subplots_adjust(left=0.05, bottom=0.02, right=1.0, top=0.98)
            slider_ax = fig.add_axes([0.05, 0.06, 0.76, 0.03])
            slider_vertical_pos = -1.05
        elif show_func in {"contourf", "contour"}:
            fig.subplots_adjust(left=0.05, bottom=0.12, right=1.05, top=0.92)
            slider_ax = fig.add_axes([0.05, 0.05, 0.895, 0.05])
            slider_vertical_pos = -0.5
        else:
            fig.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)
            slider_ax = fig.add_axes([0.05, 0.05, 0.895, 0.05])
            slider_vertical_pos = -0.5
        
        slider_ax.set_in_layout(False)
        init_val = int(0.5 * (nt - 1))
        init_percent = 0.0 if nt <= 1 else 100.0 * init_val / (nt - 1)
        frame_slider = Slider(slider_ax, "", 0, nt - 1, valinit=init_val, valstep=1)
        frame_slider.valtext.set_visible(False)
        label_text = slider_ax.text(
            0.5,
            slider_vertical_pos,
            f"Time Axis {init_percent:.0f}%",
            transform=slider_ax.transAxes,
            ha="center",
            va="center",
            fontsize=18,
            fontfamily="Times New Roman",
            color="#333333",
        )

        current_val = init_val

        def _on_slider_change(val):
            nonlocal current_val
            idx = int(val)
            current_val = idx
            update_frame(idx)
            percent = 0.0 if nt <= 1 else 100.0 * idx / (nt - 1)
            label_text.set_text(f"Time Axis {percent:.0f}%")
            fig.canvas.draw_idle()

        frame_slider.on_changed(_on_slider_change)

        def _on_key(event):
            nonlocal current_val
            if event.key in ("left", "a", "A"):
                new_val = max(0, current_val - 1)
            elif event.key in ("right", "d", "D"):
                new_val = min(nt - 1, current_val + 1)
            else:
                return
            if new_val == current_val:
                return
            frame_slider.set_val(new_val)

        fig.canvas.mpl_connect("key_press_event", _on_key)
        fig._frame_slider = frame_slider
        update_frame(init_val)
    else:
        update_frame(0)

    fig.canvas.draw_idle()

    return fig


def save_evolution_2d(
    rho: np.ndarray,
    output_path: str | Path,
    show_func: Optional[str] = None,
    fig_name: str = "Density evolution",
    barrier_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    cmap: str | None = None,
    fps: int | None = None,
    dpi: int = 120,
    codec: str | None = None,
    bitrate: int | None = None,
) -> str:
    """
    Export 2D density evolution to an animation file (e.g. .gif, .mp4).

    Notes
    -----
    - This does not change the interactive slider behavior; it is a separate export-only path.
    - When `fps` is None, the export is time-normalized to 5s total with a 0.5s pause on the first
      and last frame.
    - GIF export requires Pillow; MP4 export requires ffmpeg available on PATH.
    """
    fig, update_frame, nt, _show_func = _prepare_evolution_2d_plot(
        rho=rho,
        show_func=show_func,
        fig_name=fig_name,
        barrier_fn=barrier_fn,
        cmap=cmap,
        show_colorbar=False,
    )

    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = out_path.suffix.lower()
    if suffix not in {".gif", ".mp4", ".webm"}:
        raise ValueError("output_path must end with one of: .gif, .mp4, .webm")
    if nt <= 0:
        raise ValueError("rho must have nt > 0 to export an animation")
    if fps is None:
        # Default export: 5s total duration, with 0.5s pause at the beginning and the end.
        # We do this by resampling the time axis to a fixed number of animation frames.
        fps = int(np.clip(round(nt / 4), 2, 60))
        if fps % 2 == 1:
            fps = min(60, fps + 1)
        hold_frames = fps // 2
        core_frames = 4 * fps
        core_indices = np.rint(np.linspace(0, nt - 1, core_frames)).astype(int)
        frame_indices = [0] * hold_frames + core_indices.tolist() + [nt - 1] * hold_frames
    else:
        if fps <= 0:
            raise ValueError("fps must be positive")
        frame_indices = list(range(nt))

    # Reduce whitespace
    try:
        fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
        axes = fig.get_axes()
        if axes:
            axes[0].set_position([0, 0, 1, 1])
    except Exception:
        pass

    anim = animation.FuncAnimation(
        fig,
        func=lambda i: update_frame(int(i)),
        frames=frame_indices,
        interval=int(1000 / fps),
        blit=False,
    )

    try:
        if suffix == ".gif":
            writer = animation.PillowWriter(fps=fps)
        else:
            if not animation.FFMpegWriter.isAvailable():
                raise RuntimeError("ffmpeg is not available; cannot export mp4/webm")
            if suffix == ".webm" and codec is None:
                codec = "libvpx-vp9"
            if suffix == ".mp4" and codec is None:
                codec = "libx264"
            writer = animation.FFMpegWriter(fps=fps, codec=codec, bitrate=bitrate)
        anim.save(out_path, writer=writer, dpi=dpi)
    finally:
        plt.close(fig)

    return str(out_path)


def _adjust_fig(fig):
    """
    Adjust figure properties for better visualization.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to adjust
    """
    font_name = 'STIXGeneral'
    
    # Apply adjustments to all axes
    for ax in fig.get_axes():
        # Set font for tick labels and title
        ax.tick_params(labelsize=18)
        ax.xaxis.label.set_size(20)
        ax.yaxis.label.set_size(20)
        ax.title.set_size(20)

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontfamily(font_name)
        
        ax.xaxis.label.set_fontfamily(font_name)
        ax.yaxis.label.set_fontfamily(font_name)
        ax.title.set_fontfamily(font_name)

        legend = ax.get_legend()
        if legend is not None:
            for text in legend.get_texts():
                text.set_fontsize(16)
                text.set_fontfamily(font_name)
        
        # Set line width
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
        

    # Adjust figure layout
    fig.patch.set_facecolor('white')
