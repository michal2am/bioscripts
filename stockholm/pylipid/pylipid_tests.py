import numpy as np
import matplotlib.pyplot as plt
from pylipid.api import LipidInteraction
from pylipid.util import check_dir


trajfile_list = ["step7_production.parts1to8.xtc"]
topfile_list = ["../6x3z_CG.pdb"]

dt_traj = None            # not needed
stride = 50

lipid = "CHOL"          # residue name in the topology.
lipid_atoms = None      # means all
cutoffs = [0.5, 0.8]    # dual-cutoff scheme for coarse-grained simulations
pdb_file_to_map = None  # need to generate pdb? or maybe charmm-gui provides?


binding_site_size = 4   # binding site should contain at least four residues.
n_top_poses = 3         # write out num. of representative bound poses for each binding site.
n_clusters = "auto"

save_dir = None         # means current
save_pose_format = "pdb"
save_pose_traj = True   # save all the bound poses in a trajectory for each binding site
save_pose_traj_format = "xtc"
timeunit = "us"
resi_offset = 0         # shift the residue index, useful for MARTINI models.
radii = None            # Martini 2.2 are defined in pylipid
fig_format = "png"
num_cpus = 8

#####################################
###### no changes needed below ######
#####################################

#### calculate lipid interactions
li = LipidInteraction(trajfile_list, topfile_list=topfile_list, cutoffs=cutoffs, lipid=lipid,
                      lipid_atoms=lipid_atoms, nprot=1, resi_offset=resi_offset,
                      timeunit=timeunit, save_dir=save_dir, stride=stride, dt_traj=dt_traj)
li.collect_residue_contacts()
li.compute_residue_duration(residue_id=None)
li.compute_residue_occupancy(residue_id=None)
li.compute_residue_lipidcount(residue_id=None)
li.show_stats_per_traj(write_log=True, print_log=True)
li.compute_residue_koff(residue_id=None, plot_data=True, fig_close=True,
                        fig_format=fig_format, num_cpus=num_cpus)
li.compute_binding_nodes(threshold=binding_site_size, print_data=False)
if len(li.node_list) == 0:
    print("*"*50)
    print("No binding site detected! Skip analysis for binding sites.")
    print("*"*50)
else:
    li.compute_site_duration(binding_site_id=None)
    li.compute_site_occupancy(binding_site_id=None)
    li.compute_site_lipidcount(binding_site_id=None)
    li.compute_site_koff(binding_site_id=None, plot_data=True, fig_close=True,
                         fig_format=fig_format, num_cpus=num_cpus)
    pose_traj, pose_rmsd_data = li.analyze_bound_poses(binding_site_id=None, pose_format=save_pose_format,
                                                       n_top_poses=n_top_poses, n_clusters=n_clusters,
                                                       fig_format=fig_format, num_cpus=num_cpus)
    # save pose trajectories
    if save_pose_traj:
        for bs_id in pose_traj.keys():
            pose_traj[bs_id].save("{}/Bound_Poses_{}/Pose_traj_BSid{}.{}".format(li.save_dir, li.lipid, bs_id,
                                                                          save_pose_traj_format))
    del pose_traj  # save memory space
    surface_area_data = li.compute_surface_area(binding_site_id=None, radii=radii, fig_format=fig_format)
    data_dir = check_dir(li.save_dir, "Dataset_{}".format(li.lipid))
    pose_rmsd_data.to_csv("{}/Pose_RMSD_data.csv".format(data_dir), index=False, header=True)
    surface_area_data.to_csv("{}/Surface_Area_data.csv".format(data_dir), index=True, header=True)
    li.write_site_info(sort_residue="Residence Time")

if pdb_file_to_map is not None:
    li.save_pymol_script(pdb_file_to_map)

#### write and save data
for item in ["Dataset", "Duration", "Occupancy", "Lipid Count", "CorrCoef"]:
    li.save_data(item=item)
for item in ["Residence Time", "Duration", "Occupancy", "Lipid Count"]:
    li.save_coordinate(item=item)
for item in ["Residence Time", "Duration", "Occupancy", "Lipid Count"]:
    li.plot(item=item, fig_close=True, fig_format=fig_format)
    li.plot_logo(item=item, fig_close=True, fig_format=fig_format)

#### plot binding site comparison.
if len(li.node_list) > 0:
    for item in ["Duration BS", "Occupancy BS"]:
        li.save_data(item=item)

        ylabel_timeunit = 'ns' if li.timeunit == "ns" else r"$\mu$s"
        ylabel_dict = {"Residence Time": "Residence Time ({})".format(ylabel_timeunit),
                       "Duration": "Duration ({})".format(ylabel_timeunit),
                       "Occupancy": "Occuoancy (100%)",
                       "Lipid Count": "Lipid Count (num.)"}

        # plot No. 1
        binding_site_IDs = np.sort(
                 [int(bs_id) for bs_id in li.dataset["Binding Site ID"].unique() if bs_id != -1])
        for item in ["Residence Time", "Duration", "Occupancy", "Lipid Count"]:
            item_values = np.array(
                      [li.dataset[li.dataset["Binding Site ID"]==bs_id]["Binding Site {}".format(item)].unique()[0]
                       for bs_id in binding_site_IDs])
            fig, ax = plt.subplots(1, 1, figsize=(len(li.node_list)*0.5, 2.6))
            ax.scatter(np.arange(len(item_values)), np.sort(item_values)[::-1], s=50, color="red")
            ax.set_xticks(np.arange(len(item_values)))
            sorted_index = np.argsort(item_values)[::-1]
            ax.set_xticklabels(binding_site_IDs[sorted_index])
            ax.set_xlabel("Binding Site ID", fontsize=12)
            ax.set_ylabel(ylabel_dict[item], fontsize=12)
            for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
                plt.setp(label, fontsize=12, weight="normal")
            plt.tight_layout()
            plt.savefig("{}/{}_{}_v_binding_site.{}".format(li.save_dir, li.lipid, "_".join(item.split()), fig_format),
                        dpi=200)
            plt.close()

        # plot No. 2
        binding_site_IDs_RMSD = np.sort([int(bs_id) for bs_id in binding_site_IDs
                                        if f"Binding Site {bs_id}" in pose_rmsd_data.columns])
        RMSD_averages = np.array(
                     [pose_rmsd_data[f"Binding Site {bs_id}"].dropna(inplace=False).mean()
                      for bs_id in binding_site_IDs_RMSD])
        fig, ax = plt.subplots(1, 1, figsize=(len(li.node_list)*0.5, 2.6))
        ax.scatter(np.arange(len(RMSD_averages)), np.sort(RMSD_averages)[::-1], s=50, color="red")
        ax.set_xticks(np.arange(len(RMSD_averages)))
        sorted_index = np.argsort(RMSD_averages)[::-1]
        ax.set_xticklabels(binding_site_IDs_RMSD[sorted_index])
        ax.set_xlabel("Binding Site ID", fontsize=12)
        ax.set_ylabel("RMSD (nm)", fontsize=12)
        for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
            plt.setp(label, fontsize=12, weight="normal")
        plt.tight_layout()
        plt.savefig("{}/{}_RMSD_v_binding_site.{}".format(li.save_dir, li.lipid, fig_format), dpi=200)
        plt.close()

        # plot No. 3
        surface_area_averages = np.array(
                       [surface_area_data["Binding Site {}".format(bs_id)].dropna(inplace=False).mean()
                        for bs_id in binding_site_IDs])
        fig, ax = plt.subplots(1, 1, figsize=(len(li.node_list)*0.5, 2.6))
        ax.scatter(np.arange(len(surface_area_averages)), np.sort(surface_area_averages)[::-1], s=50, color="red")
        ax.set_xticks(np.arange(len(surface_area_averages)))
        sorted_index = np.argsort(surface_area_averages)[::-1]
        ax.set_xticklabels(binding_site_IDs[sorted_index])
        ax.set_xlabel("Binding Site ID", fontsize=12)
        ax.set_ylabel(r"Surface Area (nm$^2$)", fontsize=12)
        for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
            plt.setp(label, fontsize=12, weight="normal")
        plt.tight_layout()
        plt.savefig("{}/{}_surface_area_v_binding_site.{}".format(li.save_dir, li.lipid, fig_format), dpi=200)
        plt.close()

        # plot No. 4
        res_time_BS = np.array(
                  [li.dataset[li.dataset["Binding Site ID"]==bs_id]["Binding Site Residence Time"].unique()[0]
                   for bs_id in binding_site_IDs_RMSD])
        fig, ax = plt.subplots(1, 1, figsize=(len(li.node_list)*0.5, 2.6))
        ax.scatter(res_time_BS, RMSD_averages, s=50, color="red")
        ax.set_xlabel(ylabel_dict["Residence Time"], fontsize=12)
        ax.set_ylabel("RMSD (nm)", fontsize=12)
        for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
            plt.setp(label, fontsize=12, weight="normal")
        plt.tight_layout()
        plt.savefig("{}/{}_Residence_Time_v_RMSD.{}".format(li.save_dir, li.lipid, fig_format), dpi=200)
        plt.close()

        # plot No. 5
        res_time_BS = np.array(
                  [li.dataset[li.dataset["Binding Site ID"]==bs_id]["Binding Site Residence Time"].unique()[0]
                   for bs_id in binding_site_IDs])
        fig, ax = plt.subplots(1, 1, figsize=(len(li.node_list)*0.5, 2.6))
        ax.scatter(res_time_BS, surface_area_averages, s=50, color="red")
        ax.set_xlabel(ylabel_dict["Residence Time"], fontsize=12)
        ax.set_ylabel(r"Surface Area (nm$^2$)", fontsize=12)
        for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
            plt.setp(label, fontsize=12, weight="normal")
        plt.tight_layout()
        plt.savefig("{}/{}_Residence_Time_v_surface_area.{}".format(li.save_dir, li.lipid, fig_format), dpi=200)
        plt.close()