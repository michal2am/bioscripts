# python 3
# script for plotting vmd average lipid markers positions
# michaladammichalowski@gmail.com
# 30.11.15 - creation
# EXAMPLE CALL: python3 mm_analysis_plotthc.py --psf step5_assembly.xplor_ext.psf --dcd step7.1_production.dcd -s "resname POPC and name P" "resname POPC and name N" "resname POPC and name C1 C2 C3" --position_files "POPC_P" "POPC_N" "POPC_C" --position_plot "POPC_profile" --labels "phosphate phosphorus" "choline sodium" "glycerol carbons"

import subprocess
import argparse
import numpy as np
import mm_lib_plots as mmplt


def tcl_trajectory(psf, dcd):
    """
    :param psf: psf file path
    :param dcd: dcd file path
    :return:
    """
    script = """
set mol [mol new {0} type psf waitfor all]
mol addfile {1} type dcd step 5 waitfor all molid $mol
""".format(psf, dcd)
    return script


def tcl_avpos_z(position_file, selection):
    """
    :param position_file: name to save atom average coordinates
    :param selection: selection of atoms
    :return:
    TODO: separate upper and bottom leaflet
    """
    script = """
set file [open {0} "w"]
set selection [atomselect top "{1}"]
set positions [measure avpos $selection]
foreach sublist $positions {{
set z_position [lindex $sublist 2]
puts $file $z_position
}}
""".format(position_file, selection)
    return script


def tcl_allpos_z(position_file, selection):
    """
    :param position_file:
    :param selection:
    :return:
    """

    # alternative, measures for WHOLE selection:
    '''
    set all_position [measure center $selection]
    set z_position [lindex $all_position 2]
    '''

    position_file += ".dat"
    script = """
set file [open {0} "w"]
set selection [atomselect top "{1}"]
set nf [molinfo top get numframes]
for {{set i 0}} {{$i < $nf}} {{incr i}} {{
$selection frame $i
set z_position [$selection get z]
puts $file $z_position
}}
""".format(position_file, selection)
    return script


def create_tcl_script(script_file, *args):
    """
    :param script_file: name to save complete script
    :param args: all scripts that shall be joined
    :return:
    """
    complete = ""
    complete_file = open(script_file, "w")
    for script in args:
        complete += script
    complete += """
exit
"""
    complete_file.write(complete)


def run_tcl(script_file):
    """
    :param script_file: complete tcl script
    :return:
    """
    subprocess.call(["vmd", "-e", script_file])


def read_thc_hist(thc_files):
    """
    :param args:
    :return:
    """
    thcs = []  # container for z positions
    for thc_file in thc_files:
        thc = open(thc_file+".dat").readlines()
        par_thc = [[float(col) for col in row.split()] for row in thc]  # first: create list rows
        thcs.append([item for sublist in par_thc for item in sublist])  # second: put rows alltogether
    return thcs


def para_thc_hist(parseds_thc):
    """
    :param parseds_thc:
    :return:
    """
    for profile, selection in zip(parseds_thc, args.labels):
        top_profile = [position for position in profile if position > 0]
        bot_profile = [position for position in profile if position < 0]
        print("Density peaks of {0}: {1:.1f} SD: {2:.1f} and {3:.1f} SD: {4:.1f}, thickness: {5:.1f} SD: {6:.3f} ".format(selection, np.mean(top_profile), np.std(top_profile), np.mean(bot_profile), np.std(bot_profile), np.mean(top_profile) - np.mean(bot_profile),  (np.std(top_profile)**2 + np.std(bot_profile)**2)**0.5))



def plot_thc_hist(parseds_thc):
    """
    :param parseds_thc:
    :return:
    """
    print(len(parseds_thc))
    print(len(args.labels))
    mmplt.plot_histogram(parseds_thc, "membrane axis coordinate [A]", "density [norm]", args.labels, args.position_plot)


parser = argparse.ArgumentParser()
parser.add_argument("--psf", help="psf file name")
parser.add_argument("--dcd", help="dcd file name")
parser.add_argument("-s", "--selections", nargs='+', help="vmd atom selections")
parser.add_argument("-l", "--labels", nargs='+', help="labels for plots")
parser.add_argument("-pf", "--position_files", nargs='+', help="names to save positions")
parser.add_argument("-pp", "--position_plot", help="names to save plots")
args = parser.parse_args()

'''
read_traj = tcl_trajectory(args.psf, args.dcd)
profiles = []
for file, selection in zip(args.position_files, args.selections):
    profiles.append(tcl_allpos_z(file, selection))
create_tcl_script("complete_script.tcl", read_traj, *profiles)
run_tcl("complete_script.tcl")
'''
parseds_thc = read_thc_hist(args.position_files)
para_thc_hist(parseds_thc)
plot_thc_hist(parseds_thc)
