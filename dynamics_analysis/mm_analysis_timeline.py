# python 3
# plots vmd timeline tml files
# michaladammichalowski@gmail.com
# 15.04.16 - creation
# EXAMPLE CALL:

import argparse as arp
import numpy as np
import mm_lib_plots as mmplt
import mm_lib_analysis as mmaln


class TimelineResidue:

    def __init__(self, resid, segname):
        """
        :param resid: VMD resid number
        :param segname: VMD segment name
        """

        self.resid = resid
        self.segname = segname

        self.properties = {}
        self.means = {}

    def __str__(self):
        """
        prints single residue info
        """
        return "Resid: {} Segname: {} Properties: {}".format(self.resid, self.segname, self.properties)

    def add_property(self, prop, ini_value):
        """
        :param prop: property to add
        :param ini_value: initial value of property
        """

        self.properties.update({prop: ini_value})

    def mean_property(self, prop):
        """
        :param prop: property to calculate mean
        :return: mean value of properties
        """
        self.means.update({prop: np.mean([line[1] for line in self.properties[prop]])})


class TimelineSegment:

    def __init__(self, timeline_files, timeline_props, timestep, frequency, name):
        """
        :param timeline_files: tml files
        :param timeline_props: tml files properties
        :param timestep: md timestep
        :param frequency: md trajectory save frequency
        """

        self.name = name
        self.timeline_files = timeline_files
        self.timeline_props = timeline_props
        self.timeframe = timestep*frequency*0.000001

        self.residues = []

        self.parse_tml()

    def __contains__(self, residue):
        """
        :param residue: check if residue present in timeline
        :return: true if present
        """

        for resi in self.residues:
            if residue.resid == resi.resid and residue.segname == resi.segname:
                return True
        return False

    def get_resi(self, residue):
        """
        :param residue: selected residue
        :return: get residue from timeline
        """

        for resi in self.residues:
            if residue.resid == resi.resid and residue.segname == resi.segname:
                return resi

    def parse_tml(self):
        """
        parse tml file
        """

        for file, prop in zip(self.timeline_files, self.timeline_props):
            with open(file) as fl:
                timeline = [line.rstrip('\n').split() for line in fl.readlines() if "#" not in line]
                timeline = [[int(line[0]), line[2], float(line[3])*self.timeframe, float(line[4])] for line in timeline]

                for line in timeline:
                    #if line[2] > 5: break  # debug
                    new_residue = TimelineResidue(line[0], line[1])
                    new_residue.add_property(prop, [line[2:4]])
                    if new_residue not in self:
                        self.residues.append(new_residue)
                    else:
                        self.get_resi(new_residue).properties[prop].append(line[2:4])

                for resi in self.residues:
                    resi.mean_property(prop)

    def get_prop(self, prop):
        """
        :param prop:
        :return:
        """

        resids = [resi.resid for resi in self.residues]
        prop_values = [[resi.means[prop] for resi in self.residues]]

        return np.array(list(zip(resids, *prop_values)))


    def print_prop(self, props):
        """
        :param props: properties to plot
        """
        resids = [resi.resid for resi in self.residues]
        prop_values = []

        for prop in props:
            prop_values.append([resi.means[prop] for resi in self.residues])

        dataseries = [np.array(list(zip(resids, *prop_values)))]

        mmplt.plot_simple_multiple_numpy(dataseries, "Residue position", prop, [self.name], self.name,
                                         sizex=2.5, sizey=1.5)

class Timelines:

    def __init__(self, files, properties, timestep, frequency, names):
        self.names = names
        self.files = files
        self.properties = properties
        self.timestep = timestep
        self.frequency = frequency

        timelines = []
        split_files = [self.files[i:i+len(self.properties)] for i in range(0, len(self.files), len(self.properties))]

        for name, segment in zip(self.names, split_files):
            timelines.append(TimelineSegment(segment, self.properties, self.timestep, self.frequency, name))

        #for timeline in timelines:
            #timeline.print_prop(["RMSF"])

        dataseries = [timeline.get_prop("RMSD") for timeline in timelines]
        mmplt.plot_simple_multiple_numpy(dataseries, "Residue position", "RMSD", self.names, "timeline_RMSD",
                                         sizex=2.5, sizey=1.5)


parser = arp.ArgumentParser()
parser.add_argument("-n", "--names", nargs="+", help="segment names")
parser.add_argument("-f", "--timeline_files", nargs="+", help="timeline files to plot")
parser.add_argument("-p", "--timeline_properties", nargs="+", help="timeline property names")
parser.add_argument("-t", "--timestep", type=float, help="timestep")
parser.add_argument("-q", "--frequency", type=int, help="dcd save frequency")
args = parser.parse_args()

#timelines = TimelineSegment(args.timeline_files, args.timeline_properties, args.timestep, args.frequency)
#timelines.print_prop(["RMSD"])

test = Timelines(args.timeline_files, args.timeline_properties, args.timestep, args.frequency, args.names)
