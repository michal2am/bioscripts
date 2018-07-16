import pandas as pd
import neo
import os


class Lab2Clamp:

    def __init__(self, header, path):
        '''
        reads all ChannelLab abf files at path and saves in one atf file
        :param header: path to header template
        :param path: path to directory
        '''

        from_files = [self.from_chl(path + abf) for abf in [name for name in os.listdir(path) if name.endswith(".abf")]]
        multi_data = pd.concat(from_files, axis=1)

        self.to_pc(header, path+"/all.atf", multi_data)

    def from_chl(self, infile):
        '''
        reads single abf file
        :param infile: abf file path
        :return: abf data in Pandas data frame
        '''

        # neo-io data structure tree decomposition
        reader = neo.io.AxonIO(infile)
        bl = reader.read()[0]
        seg = bl.segments[0]
        sig = seg.analogsignals[0]

        # data parsing
        times = sig.times.rescale('ms').magnitude
        values = sig.rescale('pA').magnitude * 100                  # fake pA for percentages
        probs = values[:, 0]                                        # dirty clear

        return pd.DataFrame(index=times, data=probs)

    def to_pc(self, header, outfile, multi_data):
        '''
        saves Pandas data frame as single atf file
        :param header: path to header template
        :param outfile: path to atf outfile
        :param multi_data: Pandas data frame with all the abfs data
        :return:
        '''

        # prepare header
        with open(header, 'r') as h:

            cols = multi_data.shape[1]

            header = h.readlines()
            header[1] = "7\t{}\n".format(1+cols)
            header[8] = header[8].rstrip() + "\t" + "\"sim #1\"\t"*cols + "\n"
            header[9] = header[9].rstrip() + "\t" + "\"Trace #1 (prob)\"\t"*cols + "\n"
            header = ''.join(header)

        # write header
        with open(outfile, 'w') as d:
            d.write(header)

        # write all data
        with open(outfile, 'a') as d:
            multi_data.to_csv(d, header=False, float_format='%.4f', sep='\t')

doTheJob = Lab2Clamp("C:/Users/mm/Desktop/header.txt", "C:/Users/mm/Desktop/test/")








