from mm_solver_ode import SolverOde


class AnalyzerODE:

    def __init__(self, models, t0, te):
        self.models = models
        self.t0, self.te = t0, te
        self.results_dynamic = self.integrate_model(False)
        self.results_steady = self.integrate_model(True)


    def integrate_model(self, equi):
        return [SolverOde(model.trmn, states_ini_concentrations, states_names, t0, te, equi) for model in self.models]
