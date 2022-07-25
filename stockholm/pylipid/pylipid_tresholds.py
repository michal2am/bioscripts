from itertools import product
import mdtraj as md
import matplotlib.pyplot as plt
from pylipid.util import get_traj_info, check_dir
import numpy as np
import matplotlib.ticker as ticker


def plot_minimum_distances(distances, times, title, fn):
    fig, ax = plt.subplots(1, 1, figsize=(3, 2.5))
    ax.plot(times, distances)
    ax.set_xlabel(r"Time ($\mu$s)")
    ax.set_ylabel("Minimum distances (nm)")
    ax.set_title(title)
    ax.set_ylim(0, 1.0)
    plt.tight_layout()
    plt.savefig(fn, dpi=200)
    plt.close()
    return


def compute_minimum_distance(traj, lipid, fig_dir, lipid_atoms=None,
                            contact_frames=10, distance_threshold=0.65):
    DIST_CONTACT_ALL = []
    traj_info, _, _ = get_traj_info(traj, lipid, lipid_atoms=lipid_atoms)
    for protein_idx in np.arange(nprot, dtype=int):
        for residue_idx, residue_atom_indices in enumerate(traj_info["protein_residue_atomid_list"][protein_idx]):
            dist_matrix = np.array([np.min(md.compute_distances(traj, np.array(list(product(residue_atom_indices, lipid_atom_indices)))), axis=1) for lipid_atom_indices in traj_info["lipid_residue_atomid_list"]])
            # plot distances
            for lipid_idx in np.arange(len(dist_matrix)):
                if sum(dist_matrix[lipid_idx] < distance_threshold) >= contact_frames:
                    DIST_CONTACT_ALL.append(dist_matrix[lipid_idx])
                    plot_minimum_distances(dist_matrix[lipid_idx], traj.time/1000000.0,
                                           "{}-{}{}".format(traj_info["residue_list"][residue_idx], lipid, lipid_idx),
                                           "{}/dist_{}_{}{}.png".format(fig_dir, traj_info["residue_list"][residue_idx],
                                                                        lipid, lipid_idx))

    return DIST_CONTACT_ALL


def plot_PDF(distance_set, num_of_bins, fn):
    fig, ax = plt.subplots(1,1)
    ax.hist(distance_set, bins=num_of_bins, density=True)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Minimum distance (nm)")
    ax.set_ylabel("Probablity Density")
    ax.set_title(lipid)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    plt.tight_layout()
    plt.savefig(fn, dpi=200)
    return


# quick and dirty code to plot histogram of contact distances
# works for single xtc, files hardcoded below

trajfile = "7qn7_CG_b3_MD1.xtc"
topfile = "../7qn7_CG_b3.pdb"
lipid = "CHOL"

lipid_atoms = None # all lipid atom/bead will be considered
nprot = 1
save_dir = "test_minimum_dist_{}".format(lipid)
fig_dir = check_dir(save_dir, "Figures_dists", print_info=False)
contact_frames = 5  # will only plot data if the contact was formed over ${contact_frames} frames.
distance_threshold = 0.65

traj = md.load(trajfile, top=topfile, stride=10)
minimum_distance_set = compute_minimum_distance(traj, lipid, fig_dir, lipid_atoms=lipid_atoms, contact_frames=10,
                                                distance_threshold=0.65)

distance_set = np.concatenate(minimum_distance_set)
num_of_bins = 1000
fig_fn = "{}/dist_distribut_{}.png".format(save_dir, lipid)
plot_PDF(distance_set, num_of_bins, fig_fn)