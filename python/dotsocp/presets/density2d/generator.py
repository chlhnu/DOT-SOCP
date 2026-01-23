import os
from dotsocp.presets.resources import _resolve_resources_path
import matplotlib.image as mpimg
import numpy as np
from scipy import ndimage


__all__ = [
    "get_predefined_examples", # List of predefined examples
    "get_example", # Predefined examples
    "get_example_from_images" # Custom examples from images
]

# =========================
# Helpers
# =========================

def _gaussian_field(xx: np.ndarray, yy: np.ndarray, cx: float, cy: float, sigma: float) -> np.ndarray:
    """Isotropic 2D Gaussian."""
    dist_sq = (xx - cx) ** 2 + (yy - cy) ** 2
    return np.exp(-dist_sq / (2.0 * sigma ** 2))


def _load_and_prepare_image(name: str, target_shape: tuple[int, int], reverse: bool = False) -> np.ndarray:
    """Load an example image from resources or project examples, then resize/pad."""
    img_path = _resolve_resources_path(name)

    img = mpimg.imread(img_path)
    img = np.asarray(img, dtype=np.float64)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=2)
    if img.max() > 1.0:
        img /= 255.0

    zoom = (
        target_shape[0] / img.shape[0],
        target_shape[1] / img.shape[1],
    )
    if not np.allclose(zoom, (1.0, 1.0)):
        img = ndimage.zoom(img, zoom, order=1)

    if reverse:
        img = 1.0 - img

    # Align with plotting convention origin="lower" by flipping rows (row 0 is at bottom).
    img = np.flipud(img)

    return np.pad(img, ((0, 1), (0, 1)), mode="constant", constant_values=0.0)


def _gene_example_from_images(
    nx: int,
    ny: int,
    image_name_0: str,
    image_name_1: str,
    reverse: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Load two bitmap silhouettes and convert them into density fields."""
    target_shape = (ny - 1, nx - 1)
    rho0 = _load_and_prepare_image(image_name_0, target_shape, reverse=reverse)
    rho1 = _load_and_prepare_image(image_name_1, target_shape, reverse=reverse)
    return rho0, rho1


def _make_grid(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (xx, yy) on [0,1]^2 with shape (nx, ny), using 'ij' indexing."""
    x_lin = np.linspace(0.0, 1.0, nx)
    y_lin = np.linspace(0.0, 1.0, ny)
    return np.meshgrid(x_lin, y_lin, indexing="ij")


def _normalize_density(rho: np.ndarray, nx: int, ny: int, lower_bound: float) -> np.ndarray:
    """Scale rho so that sum ≈ nx*ny, with a simple lower_bound shift."""
    total = np.sum(rho)
    if total <= 0.0:
        raise ValueError("Density has non-positive total mass; cannot normalize.")
    scale = (nx * ny) / total
    rho = (scale * rho + lower_bound) / (1.0 + lower_bound)
    return rho


# =========================
# Example generators
# =========================

def gene_example1(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.1: Single Gaussian to single Gaussian."""
    mu1 = 0.25
    mu2 = 0.75
    sigma = np.sqrt(0.05)

    xx, yy = _make_grid(nx, ny)

    rho0 = _gaussian_field(xx, yy, mu1, mu2, sigma)
    rho1 = _gaussian_field(xx, yy, mu2, mu1, sigma)

    return rho0, rho1


def gene_example2(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.2: Single Gaussian to four Gaussian mixture."""
    mu1 = 0.25
    mu2 = 1.0 - mu1
    sigma1 = 0.1
    sigma2 = 0.05

    xx, yy = _make_grid(nx, ny)

    rho0 = _gaussian_field(xx, yy, mu1, mu1, sigma1)
    rho1 = (
        _gaussian_field(xx, yy, mu1, mu1, sigma2)
        + _gaussian_field(xx, yy, mu1, mu2, sigma2)
        + _gaussian_field(xx, yy, mu2, mu1, sigma2)
        + _gaussian_field(xx, yy, mu2, mu2, sigma2)
    )
    return rho0, rho1


def gene_example3(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.3: Laplacian to four Gaussian mixture."""
    a1, a2 = 3.0, 5.0
    mu1 = 0.25
    mu2 = 1.0 - mu1
    sigma = 0.05

    xx, yy = _make_grid(nx, ny)

    laplacian = np.exp(np.exp(-a1 * np.abs(xx - mu1) - a2 * np.abs(yy - mu1)))
    rho0 = laplacian
    rho1 = (
        _gaussian_field(xx, yy, mu1, mu1, sigma)
        + _gaussian_field(xx, yy, mu1, mu2, sigma)
        + _gaussian_field(xx, yy, mu2, mu1, sigma)
        + _gaussian_field(xx, yy, mu2, mu2, sigma)
    )
    return rho0, rho1


def gene_example4(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.4: Quartic polynomial to four Gaussian mixture."""
    xx, yy = _make_grid(nx, ny)

    center = 0.5
    rho0 = (xx - center) ** 4 + (yy - center) ** 4

    mu1 = 0.25
    mu2 = 1.0 - mu1
    sigma = 0.05
    rho1 = (
        _gaussian_field(xx, yy, mu1, mu1, sigma)
        + _gaussian_field(xx, yy, mu1, mu2, sigma)
        + _gaussian_field(xx, yy, mu2, mu1, sigma)
        + _gaussian_field(xx, yy, mu2, mu2, sigma)
    )
    return rho0, rho1


def gene_example5(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.5: Image morphing from centaur to man"""
    return _gene_example_from_images(nx, ny, "centaur.bmp", "man.bmp", reverse=True)


def _parse_stitch_env(env_value: str | None, default: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    if not env_value:
        return default
    try:
        vals = [int(x.strip()) for x in env_value.split(",")]
        if len(vals) != 4:
            raise ValueError
        return tuple(vals)  # type: ignore[return-value]
    except Exception:
        return default

def gene_example6(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.6: 4-stitch of DOTmark images."""
    img_type = os.getenv("DOTSOCP_EXAMPLE6_IMG_TYPE", "classicimages")
    img_type = (img_type or "ClassicImages").lower()
    if img_type not in {"classicimages", "shapes"}:
        raise ValueError("img_type must be 'ClassicImages' or 'Shapes'.")

    stitch1 = _parse_stitch_env(os.getenv("DOTSOCP_EXAMPLE6_STITCH1"), (1, 2, 3, 4))
    stitch2 = _parse_stitch_env(os.getenv("DOTSOCP_EXAMPLE6_STITCH2"), (5, 6, 7, 8))
    if len(stitch1) != 4:
        raise ValueError("stitch1_indices must contain exactly 4 entries.")
    if len(stitch2) != 4:
        raise ValueError("stitch2_indices must contain exactly 4 entries.")

    type_dir = "ClassicImages" if img_type == "classicimages" else "Shapes"

    def _load_and_resize(path: str, target_shape: tuple[int, int]) -> np.ndarray:
        img_path = _resolve_resources_path(path)
        img = mpimg.imread(img_path)
        img = np.asarray(img, dtype=np.float64)
        if img.ndim == 3:
            img = img[..., :3].mean(axis=2)
        if img.max() > 1.0:
            img /= 255.0

        zoom = (
            target_shape[0] / img.shape[0],
            target_shape[1] / img.shape[1],
        )
        if not np.allclose(zoom, (1.0, 1.0)):
            img = ndimage.zoom(img, zoom, order=1)

        # Align with plotting convention origin="lower" by flipping rows (row 0 is at bottom).
        return np.flipud(img)

    def _gene_big_image(paths: tuple[str, str, str, str], nx_total: int, ny_total: int) -> np.ndarray:
        sub_x = max(nx_total // 2, 1)
        sub_y = max(ny_total // 2, 1)
        target_shape = (sub_y, sub_x)  # (rows, cols)
        imgs = [
            _load_and_resize(path, target_shape)
            for path in paths
        ]
        rho = np.block([
            [imgs[0], imgs[1]],
            [imgs[2], imgs[3]],
        ])

        if rho.shape != (ny_total, nx_total):
            zoom = (
                ny_total / rho.shape[0],
                nx_total / rho.shape[1],
            )
            rho = ndimage.zoom(rho, zoom, order=1)

        return rho

    stitch1_paths = tuple(f"DOTmark/{type_dir}/{idx}.png" for idx in stitch1)
    stitch2_paths = tuple(f"DOTmark/{type_dir}/{idx}.png" for idx in stitch2)

    rho0 = _gene_big_image(stitch1_paths, nx, ny)
    rho1 = _gene_big_image(stitch2_paths, nx, ny)

    return rho0, rho1


def gene_example7(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.7: Gaussian to Dirac delta mixture."""
    # -- Random Dirac locations
    # num = 30
    # origin = [0.5, 0.5]
    # radius = 0.48 * np.random.rand(num)
    # theta = 2.0 * np.pi * np.random.rand(num)
    # dirac_x = origin[0] + radius * np.cos(theta)
    # dirac_y = origin[1] + radius * np.sin(theta)

    # -- A fixed sample for reproducibility
    dirac_x = np.array([
        0.8323, 0.5339, 0.4031, 0.6536, 0.8200, 0.4918, 0.5108, 0.6082,
        0.4633, 0.1500, 0.7227, 0.4967, 0.5318, 0.6625, 0.4309, 0.1076,
        0.3052, 0.4113, 0.4955, 0.4485, 0.5031, 0.7529, 0.4723, 0.3668,
        0.4848, 0.5474, 0.3867, 0.3192, 0.0676, 0.2382,
    ])
    dirac_y = np.array([
        0.4477, 0.6033, 0.4264, 0.5378, 0.8026, 0.7535, 0.3472, 0.2628,
        0.4023, 0.4676, 0.4535, 0.5105, 0.5903, 0.6705, 0.5134, 0.4471,
        0.6960, 0.5068, 0.5040, 0.5468, 0.2641, 0.1783, 0.2195, 0.3484,
        0.5056, 0.3925, 0.4511, 0.2659, 0.4157, 0.8016,
    ])

    xx, yy = _make_grid(nx, ny)

    rho0 = _gaussian_field(xx, yy, 0.5, 0.5, 0.1)

    rho1 = np.zeros((nx, ny), dtype=np.float64)
    hx = 1.0 / (nx - 1)
    hy = 1.0 / (ny - 1)
    x_idx = np.clip(np.round(dirac_x / hx).astype(int), 0, nx - 1)
    y_idx = np.clip(np.round(dirac_y / hy).astype(int), 0, ny - 1)
    rho1[x_idx, y_idx] = 1.0

    return rho0, rho1


def gene_example8(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.8: Highly concentrated Gaussian transportation.

    Used with the barrier of heart (see ./presets/weight2d/barrier_generator.py).
    """
    center1 = np.array([0.7, 0.3], dtype=np.float64)
    center2 = np.array([0.345, 0.625], dtype=np.float64)
    r1 = 0.09
    r2 = 0.09

    sigma1 = r1 / 3.0
    sigma2 = r2 / 3.0

    xx, yy = _make_grid(nx, ny)

    rho0 = _gaussian_field(xx, yy, center1[0], center1[1], sigma1)
    mask0 = (xx - center1[0]) ** 2 + (yy - center1[1]) ** 2 > r1 ** 2
    rho0[mask0] = 0.0

    rho1 = _gaussian_field(xx, yy, center2[0], center2[1], sigma2)
    mask1 = (xx - center2[0]) ** 2 + (yy - center2[1]) ** 2 > r2 ** 2
    rho1[mask1] = 0.0

    return rho0, rho1


def gene_example9(nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """Example 5.9: Highly concentrated Gaussian from corner to corner.
    
    Used with the barrier of maze (see ./presets/weight2d/barrier_generator.py).
    """
    center1 = np.array([0.925, 0.075], dtype=np.float64)
    center2 = np.array([0.075, 0.925], dtype=np.float64)
    r1 = 0.09
    r2 = 0.09

    sigma1 = r1 / 3.0
    sigma2 = r2 / 3.0

    xx, yy = _make_grid(nx, ny)

    rho0 = _gaussian_field(xx, yy, center1[0], center1[1], sigma1)
    mask0 = (xx - center1[0]) ** 2 + (yy - center1[1]) ** 2 > r1 ** 2
    rho0[mask0] = 0.0

    rho1 = _gaussian_field(xx, yy, center2[0], center2[1], sigma2)
    mask1 = (xx - center2[0]) ** 2 + (yy - center2[1]) ** 2 > r2 ** 2
    rho1[mask1] = 0.0

    return rho0, rho1

# =========================
# Public API
# =========================

_EXAMPLE_REGISTRY: dict[str, callable] = {
    "example1": gene_example1,
    "example2": gene_example2,
    "example3": gene_example3,
    "example4": gene_example4,
    "example5": gene_example5,
    "example6": gene_example6,
    "example7": gene_example7,
    "example8": gene_example8,
    "example9": gene_example9,
}


def get_example(problem_name: str, nx: int, ny: int, lower_bound: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (rho0, rho1) for the given example name, scaled to mass nx*ny."""
    func = _EXAMPLE_REGISTRY.get(problem_name.lower().replace('_', '-'))
    if func is None:
        # print(f"Warning: Problem '{problem_name}' not implemented, using example1.")
        # func = gene_example1
        raise ValueError(f"Problem '{problem_name}' not supported yet.")

    rho0, rho1 = func(nx, ny)

    rho0 = _normalize_density(rho0, nx, ny, lower_bound)
    rho1 = _normalize_density(rho1, nx, ny, lower_bound)

    return np.asfortranarray(rho0), np.asfortranarray(rho1)


def get_predefined_examples():
    return list(_EXAMPLE_REGISTRY.keys())


def get_example_from_images(path_to_img_0: str, path_to_img_1: str, nx: int, ny: int, lower_bound: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return (rho0, rho1) from two bitmap images, scaled to mass nx*ny."""
    rho0, rho1 = _gene_example_from_images(nx, ny, path_to_img_0, path_to_img_1)
    rho0 = _normalize_density(rho0, nx, ny, lower_bound)
    rho1 = _normalize_density(rho1, nx, ny, lower_bound)
    return np.asfortranarray(rho0), np.asfortranarray(rho1)

