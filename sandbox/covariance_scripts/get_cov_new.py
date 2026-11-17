import numpy as np
import endf
import matplotlib.pyplot as plt
import sys

import matplotlib.cm as cm
# Add the folder containing your module
sys.path.append("/Users/fosterd/Cov_getter")

# Now you can import it
import endf_cov


def plot_covariance_matrix(Erow, Ecol, C_abs, mode="absolute", info=None):
    """
    Plot a covariance matrix safely for square or rectangular matrices.

    Parameters
    ----------
    Erow : array_like
        Energies corresponding to rows (X-axis).
    Ecol : array_like
        Energies corresponding to columns (Y-axis).
    C_abs : 2D array
        Absolute covariance matrix.
    mode : str
        "absolute" for covariance, "correlation" for correlation matrix.
    info : dict
        Optional info dict to include in title (MT, MT1, LB).
    """

    # --- Determine plotting matrix ---
    if mode == "absolute":
        M = C_abs
        label = r"Covariance [$barn^2$]"
    elif mode == "correlation":
        # Handle rectangular matrices safely
        nrow, ncol = C_abs.shape

        # Compute row/column standard deviations
        sigma_row = np.sqrt(np.diag(C_abs @ C_abs.T))  # length nrow
        sigma_col = np.sqrt(np.diag(C_abs.T @ C_abs))  # length ncol

        # Avoid division by zero
        sigma_row[sigma_row == 0] = 1.0
        sigma_col[sigma_col == 0] = 1.0

        # Correlation matrix
        M = C_abs / np.outer(sigma_row, sigma_col)

        # Clip slightly outside [-1,1] due to numerical errors
        M = np.clip(M, -1.0, 1.0)
        label = "Correlation coefficient"
    else:
        raise ValueError("Unknown mode. Use 'absolute' or 'correlation'.")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(6, 5))
    pcm = ax.pcolormesh(Ecol, Erow, M, cmap="RdBu_r" if mode=="correlation" else "viridis",
                        shading="auto")
    cbar = plt.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)

    # Log-log axes
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Energy (eV)")

    # Title
    if info is None:
        ax.set_title(f"Covariance Matrix ({mode})")
    else:
        MT = info.get("MT", "")
        MT1 = info.get("MT1", "")
        LB = info.get("LB", "")
        ax.set_title(f"Covariance Matrix ({mode}), MT {MT}, MT1 {MT1}, LB {LB}")

    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    plt.show()



def code():
    Li6 = endf.Material('/Users/fosterd/nuclear_data/ENDF-B-VIII.0_neutrons/n-092_U_235.endf')
    Li6 = endf.Material('/Users/fosterd/nuclear_data/ENDF-B-VIII.0_neutrons/n-003_Li_006.endf')
    # # This is good -- leave
    for par in Li6.sections:
        if par[0] == 33:
            if Li6.section_data[par[0], par[1]]['NL'] != 0:
                # print(Li6.section_data[par[0], par[1]])
                NL = Li6.section_data[par[0], par[1]]['NL']
                for i in range(0, NL):
                    # print(Li6.section_data[par[0], par[1]]['subsections'][i])
                    MT = par[1]
                    if MT in range(850, 875):
                        print("No")
                    else:
                        # temp_driver_and_plot_MF33(par[1], Li6)
                        results = endf_cov.unpack_covariance_ni(Li6, par[1])

                        if results is not None:
                            for result in results:
                                Erow = result["Erow"]
                                Ecol = result["Ecol"]
                                C_abs = result["C_abs"]
                                info = {"MT": result["MT"], "MT1": result["MT1"], "LB": result["LB"]}
                                # Plot different views
                                plot_covariance_matrix(Erow, Ecol, C_abs, mode="absolute", info=info)
                                plot_covariance_matrix(Erow, Ecol, C_abs, mode="correlation" ,info = info)

def perturb_cross_section_KL(cross_section_energies, cross_section_values,
                             eigenvalues, eigenvectors, n_modes=None,
                             random_seed=None):
    """
    Perturbs a cross section using the KL expansion.

    Parameters:
    -----------
    cross_section_energies : ndarray
        Energy grid of the cross section (length M).
    cross_section_values : ndarray
        Original cross section values (length M).
    eigenvalues : ndarray
        Eigenvalues from covariance matrix (length N).
    eigenvectors : ndarray
        Eigenvectors interpolated onto the same energy grid (shape N x M).
    n_modes : int or None
        Number of KL modes to use (≤ N). If None, use all.
    random_seed : int or None
        For reproducibility.

    Returns:
    --------
    perturbed_cross_section : ndarray
        New cross section with KL perturbation applied (length M).
    """

    if random_seed is not None:
        np.random.seed(random_seed)

    n_total_modes = len(eigenvalues)
    M = len(cross_section_energies)

    if eigenvectors.shape != (n_total_modes, M):
        raise ValueError(f"Eigenvector array shape {eigenvectors.shape} "
                         f"does not match (n_modes, len(energies)) = ({n_total_modes}, {M})")

    if n_modes is None:
        n_modes = n_total_modes
    else:
        n_modes = min(n_modes, n_total_modes)

    # Sample standard normal random variables
    xi = np.random.randn(n_modes)

    # Clean eigenvalues
    eigenvalues_clipped = np.clip(eigenvalues[:n_modes], a_min=0.0, a_max=None)

    # KL expansion: mu(E) + sum sqrt(λ_i) * ξ_i * v_i(E)
    perturbation = np.sum(
    np.sqrt(eigenvalues_clipped[:, np.newaxis]) *xi[:, np.newaxis] * eigenvectors[:n_modes, :],
        axis=0)

    return cross_section_values + perturbation


def plot_multiple_perturbed_cross_sections_with_histogram(xs_energy,
                                                          xs_values,
                                                          eigenvalues,
                                                          eigenvectors,
                                                          n_modes,
                                                          n_perturbations,
                                                          alpha=0.5,
                                                          MT=None,
                                                          MT1=None):
    """
    Plot the original cross section and multiple KL-perturbed versions on the
    same plot, with different colors for each perturbation and a histogram of
    the spread.

    Parameters:
    - xs_energy: ndarray of energy points
    - xs_values: original cross section values
    - eigenvalues: eigenvalues of the covariance matrix
    - eigenvectors: eigenvectors (shape N x M), should already be interpolated
    - n_modes: number of KL modes to use (or None to use all)
    - n_perturbations: number of samples to draw
    - alpha: transparency level for perturbed curves
    """
    # Generate a color map to assign unique colors for each perturbation
    cmap = cm.viridis(np.linspace(0, 1, n_perturbations))

    # Store the perturbed cross sections to plot the histogram later
    perturbed_cross_sections = []

    # Set up the plot
    plt.figure(figsize=(10, 6))

    for i in range(n_perturbations):
        # Generate a perturbed cross section using KL expansion
        perturbed_xs = perturb_cross_section_KL(
            xs_energy, xs_values,
            eigenvalues=eigenvalues,
            eigenvectors=eigenvectors,
            n_modes=n_modes,
            random_seed=None  # use different random sample each time
        )

        # Plot each perturbed cross section with a unique color
        plt.plot(xs_energy, perturbed_xs, color=cmap[i],
                 alpha=alpha, linewidth=1)

        # Save the perturbed cross section for histogram
        perturbed_cross_sections.append(perturbed_xs)

    # Plot the histogram of the spread of the perturbed cross sections
    perturbed_cross_sections = np.array(perturbed_cross_sections)
    spread = np.std(perturbed_cross_sections, axis=0)

    # Plot the standard deviation (spread) as a shaded region
    plt.fill_between(xs_energy, xs_values - spread, xs_values + spread,
                     color="gray", alpha=0.3, label="Spread (Std Dev)")

    # Plot the original cross section
    plt.plot(xs_energy, xs_values, color="black", label="Original XS",
             linewidth=2)

    # Customize the plot
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Cross Section (barn)")
    plt.title(f"Original and {n_perturbations} KL-Perturbed Cross Sections with Spread. MT ="+ str(MT) + " MT1 ="+ str(MT1))
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Save the plot and show it
    plt.tight_layout()
    plt.savefig("kl_perturbed_cross_sections_with_histogram.png",
                dpi=300, bbox_inches='tight')
    plt.show()


def plot_cross_section_with_eigenvectors(xs_energy, xs_values, eigenvectors_interp, n):
    """
    Plot the cross section and the first n eigenvectors.

    Parameters:
    - xs_energy: numpy array of energy values.
    - xs_values: numpy array of cross-section values.
    - eigenvectors_interp: 2D numpy array (energy bins × eigenvector index).
    - n: number of eigenvectors to plot.
    """
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(12, 6))

    # Plot the cross section
    ax1.plot(xs_energy, xs_values, label='Cross Section', color='black', linewidth=2)
    ax1.set_xlabel("Energy (eV)")
    ax1.set_ylabel("Cross Section (barn)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Create a second y-axis for eigenvectors
    ax2 = ax1.twinx()
    for i in range(min(n, len(eigenvectors_interp))):
        ax2.plot(xs_energy, eigenvectors_interp[i], label=f"Eigenvector {i+1}", linestyle='--')
    ax2.set_ylabel("Eigenvector amplitude")
    ax2.set_yscale("linear")

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    plt.title(f"Cross Section and First {n} Eigenvectors")
    plt.tight_layout()
    plt.savefig("xsec_eigenvectors.png", dpi=300, bbox_inches='tight')
    plt.show()

def code1():
    MT = 105
    No_perturbations = 500
    Li6 = endf.Material('/Users/fosterd/nuclear_data/ENDF-B-VIII.0_neutrons/n-003_Li_006.endf')
    # Li6 = endf.Material('/Users/fosterd/nuclear_data/ENDF-B-VIII.0_neutrons/n-009_F_019.endf')

    results = endf_cov.unpack_covariance_ni(Li6, MT)
    # -- plot

    for result in results:
        Erow = result["Erow"]
        Ecol = result["Ecol"]
        C_abs =result["C_abs"]
        MT1 = result["MT1"]
        if MT1 == MT:
            info = {"MT": result["MT"], "MT1": result["MT1"], "LB": result["LB"]}
            # Plot different views
            plot_covariance_matrix(Erow, Ecol, C_abs, mode="absolute", info=info)
            plot_covariance_matrix(Erow, Ecol, C_abs, mode="correlation" ,info = info)
            # -- end plot

            eigenvalues, eigenvectors = np.linalg.eigh(C_abs)
            print(len(Erow), len(eigenvectors))

            xs_section = Li6.section_data[3, MT]['sigma']
            # E_xs, sigma_xs = xs_section.x, xs_section.y
            eigenvectors_interp = []

            for i in range(eigenvectors.shape[1]):  # loop over eigenvectors
                interp_vec = np.interp(xs_section.x, Erow[:-1], eigenvectors[:, i])

                norm = np.linalg.norm(interp_vec)
                if norm == 0:
                    print("Zero norm as eigen value = ", eigenvalues[i])
                    print("Eigenvector = ", eigenvectors[i])
                else:
                    interp_vec = interp_vec/norm

                eigenvectors_interp.append(interp_vec)

            plot_cross_section_with_eigenvectors(xs_section.x, xs_section.y,
                                                 eigenvectors_interp, n=8)

            # shape: (len(xs_energy), n_eigenvectors)
            eigenvectors_interp = np.array(eigenvectors_interp).T

            plot_multiple_perturbed_cross_sections_with_histogram(
                xs_section.x, xs_section.y,
                eigenvalues=eigenvalues,
                eigenvectors=eigenvectors_interp.T,  # shape: N x M
                n_modes=None,
                n_perturbations=No_perturbations,
                alpha=0.3,
                MT=MT,
                MT1=MT1)


if __name__ == "__main__":
    code()
