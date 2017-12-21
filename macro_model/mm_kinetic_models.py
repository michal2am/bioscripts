import logging as log
import numpy as np
from mm_kinetic import Kinetic
from mm_kinetic import Stimulus

# model definition #

class ModelBuilder:
    """
    wrapper for mm_kinetic, creates list of models of given type with selected stimulus
    """
    def __init__(self, model, agonist_concentrations, stimulus):
        """

        :param model: state/rate set selection 'jwm', 'kisiel'
        :param concentrations: list of agonist concentrations
        :param new_stimulus: stimulus type
        :return:
        """
        self.model = model
        self.stimulus = stimulus
        self.agonist_concentrations = agonist_concentrations

        self.models = []
        self.stimuli = []

        self.build_models()

        self.states = self.models[0].states
        self.states_names = self.models[0].states_names
        self.states_categories = self.models[0].states_belongs.keys()
        self.states_belongs = self.models[0].states_belongs
        self.states_number = self.models[0].states_number
        self.states_ini_concentrations = self.models[0].states_ini_concentrations

    def build_models(self):
        """
        :return:
        """
        for conc in self.agonist_concentrations:

            if self.stimulus == 'pair':
                new_stimulus = Stimulus.pair_square([0., 500.], [600., 1100.], 10.0)
                # suspend = [[5., 10.], [15., 20.]]
            if self.stimulus == 'single':
                new_stimulus = Stimulus.square(500., 1500., conc)
                # suspend = [[2.5], [10.]]

            self.stimuli.append(new_stimulus)

            log.info("### Adding rates:")

            if self.model == 'kisiel':

                r_kon   = Kinetic.State.Rate('kon',     2 * 1e7 * 1e-6,     new_stimulus)
                r_2kon  = Kinetic.State.Rate('2kon',    2 * 2 * 1e7 * 1e-6, new_stimulus)
                r_koff  = Kinetic.State.Rate('koff',    1116*1e-3)
                r_2koff = Kinetic.State.Rate('2koff',   2*1116*1e-3)
                r_b1    = Kinetic.State.Rate('b1',      150*1e-3)
                r_b2    = Kinetic.State.Rate('b2',      18000*1e-3)
                r_b3    = Kinetic.State.Rate('b3',      35*1e-3)
                r_b4    = Kinetic.State.Rate('b4',      200*1e-3)
                r_b5    = Kinetic.State.Rate('b5',      100*1e-3)
                r_a1    = Kinetic.State.Rate('a1',      20000*1e-3)
                r_a2    = Kinetic.State.Rate('a2',      800*1e-3)
                r_a3    = Kinetic.State.Rate('a3',      3333*1e-3)
                r_a4    = Kinetic.State.Rate('a4',      17500*1e-3)
                r_a5    = Kinetic.State.Rate('a5',      3000*1e-3)
                r_d1    = Kinetic.State.Rate('d1',      310*1e-3)
                r_d2    = Kinetic.State.Rate('d2',      800*1e-3)
                r_d4    = Kinetic.State.Rate('d4',      1000*1e-3)
                r_r1    = Kinetic.State.Rate('r1',      5*1e-3)
                r_r2    = Kinetic.State.Rate('r2',      400*1e-3)
                r_r4    = Kinetic.State.Rate('r4',      10*1e-3)
                r_y2    = Kinetic.State.Rate('y2',      4400*1e-3)
                r_g2    = Kinetic.State.Rate('g2',      4300*1e-3)
                r_0     = Kinetic.State.Rate('block',   0.00*1e-3)

            if self.model == 'jwm':

                r_kon   = Kinetic.State.Rate('kon',     2.00, new_stimulus)
                r_2kon  = Kinetic.State.Rate('2kon',    4.00, new_stimulus)
                r_koff  = Kinetic.State.Rate('koff',    4.00)
                r_2koff = Kinetic.State.Rate('2koff',   8.00)
                r_d     = Kinetic.State.Rate('d',       1.00)
                r_r     = Kinetic.State.Rate('r',       0.70)
                r_b     = Kinetic.State.Rate('b',       5.00)
                r_a     = Kinetic.State.Rate('a',       3.00)
                r_0     = Kinetic.State.Rate('block',   0.00)

            if self.model == 'fjwm':
                r_kon   = Kinetic.State.Rate('kon',     9.9, new_stimulus)
                r_2kon  = Kinetic.State.Rate('2kon',    19.8, new_stimulus)
                r_koff  = Kinetic.State.Rate('koff',    1.16)
                r_2koff = Kinetic.State.Rate('2koff',   2.32)
                r_d     = Kinetic.State.Rate('d',       23.8)
                r_r     = Kinetic.State.Rate('r',       0.12)
                r_b     = Kinetic.State.Rate('b',       16.5)
                r_a     = Kinetic.State.Rate('a',       1.69)
                r_y     = Kinetic.State.Rate('y',       4.46)
                r_g     = Kinetic.State.Rate('g',       4.03)
                r_mf    = Kinetic.State.Rate('mf',      0.00)
                r_mb    = Kinetic.State.Rate('mb',      0.00)
                r_sf    = Kinetic.State.Rate('sf',      0.001)
                r_sb    = Kinetic.State.Rate('sb',      0.03)
                r_0     = Kinetic.State.Rate('block',   0.00)

            if self.model == 'STOP_fjwm':
                r_kon   = Kinetic.State.Rate('kon',     0.10, new_stimulus)
                r_2kon  = Kinetic.State.Rate('2kon',    0.20, new_stimulus)
                r_koff  = Kinetic.State.Rate('koff',    1.16)
                r_2koff = Kinetic.State.Rate('2koff',   2.32)
                r_d     = Kinetic.State.Rate('d',       28.8)
                r_r     = Kinetic.State.Rate('r',       0.21)
                r_b     = Kinetic.State.Rate('b',       16.5)
                r_a     = Kinetic.State.Rate('a',       1.69)
                r_y     = Kinetic.State.Rate('y',       195)
                r_g     = Kinetic.State.Rate('g',       0.27)
                r_mf    = Kinetic.State.Rate('mf',      0.00)
                r_mb    = Kinetic.State.Rate('mb',      0.00)
                r_0     = Kinetic.State.Rate('block',   0.00)

                '''
                if self.model == 'kisiel_17':
                    r_kon = Kinetic.State.Rate('kon', 43.2, new_stimulus)
                    r_2kon = Kinetic.State.Rate('2kon', 86.4, new_stimulus)
                    r_koff = Kinetic.State.Rate('koff', 1.82)
                    r_2koff = Kinetic.State.Rate('2koff', 3.64)
                    r_a0 = Kinetic.State.Rate('a0 ', 23.8)
                    r_r = Kinetic.State.Rate('r', 0.12)
                    r_b = Kinetic.State.Rate('b', 16.5)
                    r_a = Kinetic.State.Rate('a', 1.69)
                    r_y = Kinetic.State.Rate('y', 4.46)
                    r_g = Kinetic.State.Rate('g', 4.03)
                    r_s1f = Kinetic.State.Rate('s1f', 0.001)
                    r_s1b = Kinetic.State.Rate('s1b', 15)
                    r_s2f = Kinetic.State.Rate('s2f', 0.5)
                    r_s2b = Kinetic.State.Rate('s2b', 4.0)
                    r_sdf = Kinetic.State.Rate('sdf', 1.0)
                    r_sdb = Kinetic.State.Rate('sdb', 0.5)
                    r_0 = Kinetic.State.Rate('block', 0.00)
                '''

            log.info("### Adding states:")

            if self.model == 'kisiel':

                #                                  A1O      A2O     A3O     A40     A50     R       A1R     A2R     A2F     A1D     A2D     A4D
                st_a1o  = Kinetic.State(0, 'A1O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_a1,   r_0,    r_0,    r_d1,   r_0,    r_0],   0, 'open')
                st_a2o  = Kinetic.State(1, 'A2O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a2,   r_0,    r_0,    r_0],   0, 'open')
                st_a3o  = Kinetic.State(2, 'A3O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a3,   r_0,    r_0],   0, 'open')
                st_a4o  = Kinetic.State(3, 'A4O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a4,   r_0,    r_0,    r_0,    r_d4],  0, 'open')
                st_a5o  = Kinetic.State(4, 'A5O', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_a5],  0, 'open')
                st_r    = Kinetic.State(5, 'R'  , [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_2kon, r_0,    r_0,    r_0,    r_0,    r_0],   1, 'unbound')
                st_a1r  = Kinetic.State(6, 'A1R', [r_a1,    r_0,    r_0,    r_0,    r_0,    r_koff, r_0,    r_kon,  r_0,    r_0,    r_0,    r_0],   0, 'single-bound')
                st_a2r  = Kinetic.State(7, 'A2R', [r_0,     r_0,    r_0,    r_b4,   r_0,    r_0,    r_2koff,r_0,    r_g2,   r_0,    r_0,    r_0],   0, 'double-bound')
                st_a2f  = Kinetic.State(8, 'A2F', [r_0,     r_b2,   r_0,    r_0,    r_0,    r_0,    r_0,    r_y2,    r_0,   r_0,    r_d2,   r_0],   0, 'flipped')
                st_a1d  = Kinetic.State(9, 'A1D', [r_r1,    r_0,    r_b3,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, 'desensitized')
                st_a2d  = Kinetic.State(10,'A2D', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_r2,   r_0,    r_0,    r_0],   0, 'desensitized')
                st_a4d  = Kinetic.State(11,'A4D', [r_0,     r_0,    r_0,    r_r4,   r_b5,   r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0],   0, 'desensitized')

            if self.model == 'jwm':

                st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0],   1, False)
                st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0],   0, False)
                st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_d,  r_b],   0, False)
                st_a2d  = Kinetic.State(3, 'A2D', [r_0,     r_0,      r_r,    r_0,  r_0],   0, False)
                st_a2o  = Kinetic.State(4, 'A2O', [r_0,     r_0,      r_a,    r_0,  r_0],   0, True)

            if self.model == 'STOP_fjwm':

                #                                  R        AR        A2R     A2F   A2D   A2O
                st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0,  r_0],   100, 'unbound')
                st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0,  r_0],   0, 'single-bound')
                st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_y,  r_mf, r_0],   0, 'double-bound')
                st_a2f  = Kinetic.State(3, 'A2F', [r_0,     r_2koff,  r_g,    r_0,  r_d,  r_b],   0, 'flipped')
                st_a2d  = Kinetic.State(4, 'A2D', [r_0,     r_0,      r_mb,   r_r,  r_0,  r_0],   0, 'desensitized')
                st_a2o  = Kinetic.State(5, 'A2O', [r_0,     r_0,      r_0,    r_a,  r_0,  r_0],   0, 'open')

            if self.model == 'fjwm':

                #                                  R        AR        A2R     A2F   A2D   A2O   S
                st_r    = Kinetic.State(0, 'R',   [r_0,     r_2kon,   r_0,    r_0,  r_0,  r_0,  r_sf,],   100, 'unbound')
                st_ar   = Kinetic.State(1, 'AR',  [r_koff,  r_0,      r_kon,  r_0,  r_0,  r_0,  r_0,],   0, 'single-bound')
                st_a2r  = Kinetic.State(2, 'A2R', [r_0,     r_2koff,  r_0,    r_y,  r_mf, r_0,  r_0,],   0, 'double-bound')
                st_a2f  = Kinetic.State(3, 'A2F', [r_0,     r_2koff,  r_g,    r_0,  r_d,  r_b,  r_0,],   0, 'flipped')
                st_a2d  = Kinetic.State(4, 'A2D', [r_0,     r_0,      r_mb,   r_r,  r_0,  r_0,  r_0,],   0, 'desensitized')
                st_a2o  = Kinetic.State(5, 'A2O', [r_0,     r_0,      r_0,    r_a,  r_0,  r_0,  r_0,],   0, 'open')
                st_so   = Kinetic.State(6,  'S',  [r_sb,    r_0,      r_0,    r_0,  r_0,  r_0,  r_0,],   0, 'open' )

            if self.model == 'sfjwm':

                r_kon   = Kinetic.State.Rate('kon',     9.9, new_stimulus)
                r_2kon  = Kinetic.State.Rate('2kon',    19.8, new_stimulus)
                r_koff  = Kinetic.State.Rate('koff',    1.16)
                r_2koff = Kinetic.State.Rate('2koff',   2.32)
                r_d     = Kinetic.State.Rate('d',       23.8)
                r_r     = Kinetic.State.Rate('r',       0.12)
                r_b     = Kinetic.State.Rate('b',       16.5)
                r_a     = Kinetic.State.Rate('a',       1.69)
                r_y     = Kinetic.State.Rate('y',       4.46)
                r_g     = Kinetic.State.Rate('g',       4.03)
                r_s1f   = Kinetic.State.Rate('s1f',     0.001)
                r_s1b   = Kinetic.State.Rate('s1b',     15)
                r_s2f   = Kinetic.State.Rate('s2f',     0.5)
                r_s2b   = Kinetic.State.Rate('s2b',     4.0)
                r_sdf   = Kinetic.State.Rate('sdf',     1.0)
                r_sdb   = Kinetic.State.Rate('sdb',     0.5)
                r_0     = Kinetic.State.Rate('block',   0.00)

                #                                  R         AR        A2R     A2F   A2D   A2O   RS1     RSD     RS2
                st_r    = Kinetic.State(0, 'R',    [r_0,     r_2kon,   r_0,    r_0,  r_0,  r_0,  r_s1f,  r_0,    r_0  ],   100, 'unbound')
                st_ar   = Kinetic.State(1, 'AR',   [r_koff,  r_0,      r_kon,  r_0,  r_0,  r_0,  r_0,    r_0,    r_0  ],   0, 'single-bound')
                st_a2r  = Kinetic.State(2, 'A2R',  [r_0,     r_2koff,  r_0,    r_y,  r_0,  r_0,  r_0,    r_0,    r_0  ],   0, 'double-bound')
                st_a2f  = Kinetic.State(3, 'A2F',  [r_0,     r_2koff,  r_g,    r_0,  r_d,  r_b,  r_0,    r_0,    r_0  ],   0, 'flipped')
                st_a2d  = Kinetic.State(4, 'A2D',  [r_0,     r_0,      r_0,    r_r,  r_0,  r_0,  r_0,    r_0,    r_0  ],   0, 'desensitized')
                st_a2o  = Kinetic.State(5, 'A2O',  [r_0,     r_0,      r_0,    r_a,  r_0,  r_0,  r_0,    r_0,    r_0  ],   0, 'open')
                st_rs1   = Kinetic.State(6, 'RS1', [r_s1b,   r_0,      r_0,    r_0,  r_0,  r_0,  r_0,    r_sdf,  r_0  ],   0, 'open' )
                st_rsd   = Kinetic.State(7, 'RSD', [r_0,     r_0,      r_0,    r_0,  r_0,  r_0,  r_sdb,  r_0,    r_s2f],   0, 'desensitized' )
                st_rs2   = Kinetic.State(8, 'RS2', [r_0,     r_0,      r_0,    r_0,  r_0,  r_0,  r_0,    r_s2b,  r_0  ],   0, 'open' )

            if self.model == 'spont18':

                r_b0  = Kinetic.State.Rate('b0',    0.50)                   # ?
                r_a0  = Kinetic.State.Rate('a0',    19.5)                   # sum-fit do short shut single (20.50)
                r_d0  = Kinetic.State.Rate('d0',    1.00)                   # jw
                r_r0  = Kinetic.State.Rate('r0',    0.50)                   # Kisiel17
                r_b0p = Kinetic.State.Rate('b0p',   0.25)                   # manual (10% stanów RO2 w singlu)
                r_a0p = Kinetic.State.Rate('a0p',   3.64)                   # fit do long shut single channel

                r_onm = Kinetic.State.Rate('bon',   12.00, new_stimulus)    # ok. 5e4 wysyca, nisko żeby nie było peak'u
                r_ofm = Kinetic.State.Rate('bof',   0.01)                  # musi być wolno (deaktywacja)
                r_bm  = Kinetic.State.Rate('bm',    6.00)                   # ?
                r_am  = Kinetic.State.Rate('am',    10.37)                   # sum-fit do short shut single (14.37)
                r_dm  = Kinetic.State.Rate('dm',    4.00)                   # jw
                r_rm  = Kinetic.State.Rate('rm',    0.25)                   # manual
                r_bmp = Kinetic.State.Rate('bmp',   0.01)                   # manual (10% stanów RMO2 w singlu)
                r_amp = Kinetic.State.Rate('amp',   2.15)                   # fit do long shut single channel

                r_0   = Kinetic.State.Rate('block',   0.00)

                #                                     R        RO1     RD      RO2     RM      RMO1    RMD     RMO2
                st_r      = Kinetic.State(0, 'R',    [r_0,     r_b0,   r_0,    r_0,    r_onm,  r_0,    r_0,    r_0,],   100,  'unbound')
                st_ro1    = Kinetic.State(1, 'RO1',  [r_a0,    r_0,    r_d0,   r_0,    r_0,    r_0,    r_0,    r_0,],   0,    'open')
                st_rd     = Kinetic.State(2, 'RD',   [r_0,     r_r0,   r_0,    r_b0p,  r_0,    r_0,    r_0,    r_0,],   0,    'unbound')
                st_ro2    = Kinetic.State(3, 'RO2',  [r_0,     r_0,    r_a0p,  r_0,    r_0,    r_0,    r_0,    r_0,],   0,    'open')
                st_rm     = Kinetic.State(4, 'RM',   [r_ofm,   r_0,    r_0,    r_0,    r_0,    r_bm,   r_0,    r_0,],   0,    'unbound')
                st_rmo1   = Kinetic.State(5, 'RMO1', [r_0,     r_0,    r_0,    r_0,    r_am,   r_0,    r_dm,   r_0,],   0,    'open')
                st_rmd    = Kinetic.State(6, 'RMD',  [r_0,     r_0,    r_0,    r_0,    r_0,    r_rm,   r_0,    r_bmp,], 0,    'unbound')
                st_rmo2   = Kinetic.State(7, 'RMO2', [r_0,     r_0,    r_0,    r_0,    r_0,    r_0,    r_amp,  r_0,],   0,    'open')

                self.states = [st_r, st_ro1, st_rd, st_ro2, st_rm, st_rmo1, st_rmd, st_rmo2]
                model_kinetic = Kinetic(self.states)

            '''
            if self.model == 'kisiel_17':

                #                                       SO1     SO2     A1O1    A2O2    A2OSH1  A2OSH2  A2O     R       SC      A1R     A1F     A1D     A2R     A2DSH   A2F     A2D1    A2D2
                st_so1    = Kinetic.State(0, 'SO1',    [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_so2    = Kinetic.State(0, 'SO2',    [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a1o1   = Kinetic.State(0, 'A1O1',   [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a1o2   = Kinetic.State(0, 'A1O2',   [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a2osh1 = Kinetic.State(0, 'A2OSH1', [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      ?,      r_0,    r_0,    r_0,],    0, 'open')
                st_a2osh2 = Kinetic.State(0, 'A2OSH2', [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,],    0, 'open')
                st_a2o    = Kinetic.State(0, 'A2O',    [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,],    0, 'open')
                st_r      = Kinetic.State(0, 'R',      [?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_sc     = Kinetic.State(0, 'SC',     [?,      ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a1r    = Kinetic.State(0, 'A1R',    [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    ?,      r_0,    ?,      r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a1f    = Kinetic.State(0, 'A1F',    [r_0,    r_0,    ?,      ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      ?,      r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a1d    = Kinetic.State(0, 'A1D',    [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a2r    = Kinetic.State(0, 'A2R',    [r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,],    0, 'open')
                st_a2dsh  = Kinetic.State(0, 'A2DSH',  [r_0,    r_0,    r_0,    r_0,    ?,      ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,],    0, 'open')
                st_a2f    = Kinetic.State(0, 'A2F',    [r_0,    r_0,    r_0,    ?,      r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,    ?,       ?, ],    0, 'open')
                st_a2d1   = Kinetic.State(0, 'A2D1',   [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,],    0, 'open')
                st_a2d2   = Kinetic.State(0, 'A2D2',   [r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    r_0,    ?,      r_0,    r_0,],    0, 'open')
            '''

            if self.model == 'kisiel':

                self.states = [st_a1o, st_a2o, st_a3o, st_a4o, st_a5o, st_r, st_a1r, st_a2r, st_a2f, st_a1d, st_a2d, st_a4d]
                model_kinetic = Kinetic(self.states)

            if self.model == 'jwm':
                self.states = [st_r, st_ar, st_a2r, st_a2d, st_a2o]
                model_kinetic = Kinetic(self.states)

            if self.model == 'fjwm':
                self.states = [st_r, st_ar, st_a2r, st_a2f, st_a2d, st_a2o, st_so]
                model_kinetic = Kinetic(self.states)

            if self.model == 'sfjwm':
                self.states = [st_r, st_ar, st_a2r, st_a2f, st_a2d, st_a2o, st_rs1, st_rsd, st_rs2]
                model_kinetic = Kinetic(self.states)

            self.models.append(model_kinetic)

