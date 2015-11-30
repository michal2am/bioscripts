# python 3
# script for plotting vmd average lipid markers positions
# michaladammichalowski@gmail.com
# 30.11.15 - creation

import subprocess
import argparse
import mm_lib_plots as mmplt


def create_tcl_trajectory(psf, dcd):
    """
    :param psf: psf file path
    :param dcd: dcd file path
    :return:
    """
    script = """
set mol [mol new {0} type psf waitfor all]
mol addfile {1} type dcd step 10 waitfor all molid $mol
""".format(psf, dcd)
    return script


def create_tcl_selection(position_file, selection):
    """
    :param position_file: name to save atom average coordinates
    :param selection: selection of atoms
    :return:
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


def read_thc(*args):
    """
    :param args:
    :return:
    """

    '''
    dist_rdf = []  # container for distances x
    dens_rdf = []  # container for densities y

    for rad_file in rad_files:
        rdf = open(rad_file).readlines()
        par_rdf = [[float(col) for col in row.split()] for row in rdf[0:-1]]  # first: create list rows
        dist_rdf.append([row[0] for row in par_rdf])  # second: select x columns from rows and add to container
        dens_rdf.append([row[1] for row in par_rdf])  # second: select y columns from rows and add to container

    return [dist_rdf, dens_rdf]
    '''
    pass

def plot_thc():
    """
    :return:
    """

    '''
    mmplt.plot_simple_multiple(parseds_rdf[0], parseds_rdf[1], "distance [A]", "density [rel]", args.labels, args.out_file)
    '''
    pass

parser = argparse.ArgumentParser()
parser.add_argument("--psf", help="psf file name")
parser.add_argument("--dcd", help="dcd file name")
parser.add_argument("-s", "--selections", nargs='+', help="vmd atom selections")
parser.add_argument("-pf", "--position_files", nargs='+', help="names to save positions")
parser.add_argument("-pp", "--position_plots", nargs='+', help="names to save plots")
args = parser.parse_args()

read_traj = create_tcl_trajectory(args.psf, args.dcd)
prof1 = create_tcl_selection("POPC_P_top.dat", "resname POPC and name P and x > 0")
prof2 = create_tcl_selection("POPC_P_bot.dat", "resname POPC and name P and x < 0")
create_tcl_script("complete_script.tcl", read_traj, prof1, prof2)
run_tcl("complete_script.tcl")
