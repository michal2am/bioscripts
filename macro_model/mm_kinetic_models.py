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
                new_stimulus = Stimulus.pair_square([0., 10.], [5., 15.], 10.0)
                # suspend = [[5., 10.], [15., 20.]]
            if self.stimulus == 'single':
                new_stimulus = Stimulus.square(0., 25., conc)
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

            if self.model == 'kisiel':

                self.states = [st_a1o, st_a2o, st_a3o, st_a4o, st_a5o, st_r, st_a1r, st_a2r, st_a2f, st_a1d, st_a2d, st_a4d]
                model_kinetic = Kinetic(self.states)

            if self.model == 'jwm':
                self.states = [st_r, st_ar, st_a2r, st_a2d, st_a2o]
                model_kinetic = Kinetic(self.states)

            self.models.append(model_kinetic)

