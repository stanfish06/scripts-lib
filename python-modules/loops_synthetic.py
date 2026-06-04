import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from anndata import AnnData

def sampling_donut_2d(
    r_inner: float,
    r_outer: float,
    n_points: int,
    gauss_noise_std: float = 0.0,
    n_radius_wiggles: int = 0,
    amplitude_radius_wiggle: float = 0.1,
    random_seed: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample points from a 2D annulus (donut) w/wo a sinusoidally wiggled radius.

    Points are drawn uniformly over the annular area between ``r_inner`` and
    ``r_outer``. The radius is then scaled by a sine wave in the angular
    coordinate, and optional isotropic Gaussian noise is added in Cartesian
    space.

    Parameters
    ----------
    r_inner : float
        Inner radius of the annulus.
    r_outer : float
        Outer radius of the annulus.
    n_points : int
        Number of points to sample.
    gauss_noise_std : float, default 0.0
        Standard deviation of the Gaussian noise added independently to the
        x and y coordinates.
    n_radius_wiggles : int, default 5
        Number of sinusoidal oscillations of the radius over one full turn.
    amplitude_radius_wiggle : float, default 0.1
        Relative amplitude of the radius oscillation, as a fraction of the
        local radius.
    random_seed : int, default 1
        Seed for the NumPy global random state.

    Returns
    -------
    coords : np.ndarray
        Array of shape ``(n_points, 2)`` with the sampled (x, y) coordinates.
    theta : np.ndarray
        Array of shape ``(n_points,)`` with the angular coordinate of each
        point.
    """

    np.random.seed(random_seed)
    """ sample points evenly in polar coordinates
    𝛉 ~ 𝒰(0, 2𝛑)
    R does not follow uniform otherwise points will be denser close to center
    CDF of R should be P(R < r) = (r^2 - r_inner^2) / (r_outer^2 - r_inner^2) (proportional to the area)
    then let P(R < r) = u ~ 𝒰(0, 1), r = sqrt(u * (r_max^2 - r_min^2) + r_inner^2)
    """
    r = np.sqrt(np.random.rand(n_points) * (r_outer**2 - r_inner**2) + r_inner**2)
    theta = 2 * np.pi * np.random.rand(n_points)
    # wiggle is proportional to r
    wiggle = (amplitude_radius_wiggle * np.sin(n_radius_wiggles * theta)) * r
    r_wiggled = r + wiggle
    x = r_wiggled * np.cos(theta) + np.random.normal(0, gauss_noise_std, n_points)
    y = r_wiggled * np.sin(theta) + np.random.normal(0, gauss_noise_std, n_points)
    coords = np.column_stack([x, y])
    return coords, theta


def computing_isometric_embedding(
    data: np.ndarray,
    target_dim: int,
    gauss_noise_std: float = 0.1,
    random_seed: int = 1,
) -> np.ndarray:
    """Embed low-dimensional data into a higher-dimensional space isometrically.

    The input is mapped through a random orthonormal projection, after which isotropic
    Gaussian noise is added in the target space.

    Parameters
    ----------
    data : np.ndarray
        Input array of shape ``(n_points, dim_data)``.
    target_dim : int
        Dimension of the embedding space. Must be at least ``dim_data``.
    gauss_noise_std : float, default 0.1
        Standard deviation of the Gaussian noise added in the target space.
    random_seed : int, default 1
        Seed for the NumPy global random state.

    Returns
    -------
    embedded : np.ndarray
        Array of shape ``(n_points, target_dim)`` with the embedded data.
    """
    np.random.seed(random_seed)
    n_points, dim_data = data.shape
    # n random vectors in the n-dim space (might be dependent but most should be independent)
    random_matrix = np.random.randn(target_dim, target_dim)
    # obtain orthogonal bases in the column space of random_matrix
    Q, _ = np.linalg.qr(random_matrix)
    # subspace to which the low-d data will be projected
    subspace_projection = Q[:, :dim_data]
    embedded_signal = data @ subspace_projection.T
    gaussian_noise = np.random.randn(n_points, target_dim) * gauss_noise_std
    embedded = embedded_signal + gaussian_noise
    return embedded


def main() -> dict[int, AnnData]:
    r_inner = 1.0
    r_outer = 2.0
    n_points = 2000
    donut_noise_std = 0.05
    n_radius_wiggles = 5
    amplitude_radius_wiggle = 0.1
    random_seed = 1
    target_dims = [10, 50, 100]
    embedding_noise_std = 0.1

    # Generate 2D donut with wiggles
    print(
        f"Generating 2D donut with {n_radius_wiggles} wiggles "
        f"(amplitude={amplitude_radius_wiggle})..."
    )
    coords_2d, theta = sampling_donut_2d(
        r_inner=r_inner,
        r_outer=r_outer,
        n_points=n_points,
        gauss_noise_std=donut_noise_std,
        n_radius_wiggles=n_radius_wiggles,
        amplitude_radius_wiggle=amplitude_radius_wiggle,
        random_seed=random_seed,
    )

    # Visualize the wiggled donut
    plt.figure(figsize=(6, 6))
    plt.scatter(
        coords_2d[:, 0],
        coords_2d[:, 1],
        c=theta,
        cmap="twilight",
        s=5,
        alpha=0.7,
    )
    plt.title(
        f"Wiggled Donut (n_wiggles={n_radius_wiggles}, "
        f"amplitude={amplitude_radius_wiggle})"
    )
    plt.axis("equal")
    plt.tight_layout()
    plt.show()

    # Create AnnData objects with isometric embeddings
    adata_dict: dict[int, AnnData] = {}
    for dim in target_dims:
        print(
            f"\nCreating {dim}D embedding "
            f"(with Gaussian noise scale={embedding_noise_std})..."
        )
        embedded = computing_isometric_embedding(
            coords_2d,
            target_dim=dim,
            gauss_noise_std=embedding_noise_std,
            random_seed=random_seed,
        )

        adata = AnnData(embedded)
        adata.obsm["X_original_2d"] = coords_2d
        adata.obs["theta"] = theta
        adata.uns["embedding_dim"] = dim
        adata.uns["noise_scale"] = embedding_noise_std
        adata.uns["n_wiggles"] = n_radius_wiggles
        adata.uns["wiggle_amplitude"] = amplitude_radius_wiggle

        # Compute dimensionality reductions
        print("  Computing PCA...")
        n_pca_comps = min(50, dim - 1, n_points - 1)  # Must be < min(n_samples, n_features)
        sc.pp.pca(adata, n_comps=n_pca_comps)

        print("  Computing neighbors...")
        sc.pp.neighbors(adata, n_neighbors=30, use_rep="X")

        print("  Computing UMAP...")
        sc.tl.umap(adata)

        print("  Computing t-SNE...")
        sc.tl.tsne(adata, n_pcs=n_pca_comps)

        print("  Computing Diffusion Map...")
        sc.pp.neighbors(adata, n_neighbors=100, method="gauss")
        sc.tl.diffmap(adata)

        adata_dict[dim] = adata

    return adata_dict

if __name__ == "__main__":
    main()
