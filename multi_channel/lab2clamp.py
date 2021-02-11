# concatenates channel lab traces into multisweep files

import pandas as pd
import numpy as np
import neo
import os


class Lab2Clamp:

    def __init__(self, header, fpath):
        '''
        reads all ChannelLab abf files at path and saves in one atf file
        :param header: path to header template
        :param path: path to directory
        '''

        files = os.listdir(fpath)
        files.sort(key=lambda s: os.path.getmtime(os.path.join(fpath, s)))

        #file_names = [name for name in sorted(os.listdir(fpath), key=os.path.getmtime) if name.endswith(".abf")]

        from_files = [self.from_chl(fpath + abf) for abf in files]
        multi_data = pd.concat(from_files, axis=1)

        self.to_pc(header, fpath+"/all.atf", multi_data)

    def from_chl(self, infile):
        '''
        reads single abf file
        :param infile: abf file path
        :return: abf data in Pandas data frame
        '''

        print("parsing file: {}".format(infile))

        # neo-io data structure tree decomposition
        reader = neo.io.AxonIO(infile)
        bl = reader.read()[0]
        seg = bl.segments[0]
        sig = seg.analogsignals[0]

        # data parsing
        times = sig.times.rescale('ms').magnitude
        values = sig.rescale('pA').magnitude * 100                  # fake pA for percentages
        probs = values[:, 0]                                        # dirty clear

        # dirty pClamp fit for x axis extension at time shifting
        #times = np.append(times, [times[1]*i + times[-1] for i in range(1, 100000)])
        #probs = np.append(probs, [0 for i in range(1, 100000)])

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


doTheJob = Lab2Clamp(r"C:\Users\mm\Desktop\header.txt", r"C:\Users\mm\Google Drive\RESEARCH\WORK\new_model_F200\kinetic\tyr_500\\")