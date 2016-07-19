import numpy as np
import dcprogs as dc

tau = 1e-4
tcritical = 5e-3

matrix = [ [-3050,        50,  3000,      0,    0],
           [2./3., -1502./3.,     0,    500,    0],
           [   15,         0, -2065,     50, 2000],
           [    0,     15000,  4000, -19000,    0],
           [    0,         0,    10,      0,  -10] ]

qmatrix = dc.likelihood.QMatrix(matrix, 2)
print(qmatrix)

bursts = [  [0.1, 0.2, 0.1],                  # 1st burst
            [0.2],                            # 2nd burst
            [0.15, 0.16, 0.18, 0.05, 0.1] ]   # 3rd burst

likelihood = dc.likelihood.Log10Likelihood(bursts, nopen=2, tau=tau, tcritical=tcritical)
print(likelihood)

result = likelihood(matrix)
print("Computation: {0}".format(result))

eG = dc.likelihood.MissedEventsG(qmatrix, tau=tau)
idealG = dc.likelihood.IdealG(qmatrix)

print("Equilibrium Occupancies\n"            \
      "=======================\n\n"          \
      "Ideal Likelihood\n"                   \
      "----------------\n\n"                 \
      "  * initial: {ideal_initial!r}\n"     \
      "  * final: {ideal_final!r}\n\n\n"     \
      "Missed-events Likelihood\n"           \
      "------------------------\n\n"         \
      "  * initial: {equi_initial!r}\n"      \
      "  * final: {equi_final!r}\n\n\n\n"    \
      "CHS Occupancies\n"                    \
      "===============\n\n"                  \
      "Missed-events Likelihood\n"           \
      "------------------------\n\n"         \
      "  * tcritical: {tcritical}\n"         \
      "  * initial: {chs_initial!r}\n"       \
      "  * final: {chs_final!r}"             \
      .format(
        ideal_initial = idealG.initial_occupancies,
        ideal_final   = idealG.final_occupancies,
        equi_initial  = eG.initial_occupancies,
        equi_final    = eG.final_occupancies,
        chs_initial   = eG.initial_CHS_occupancies(tcritical),
        chs_final     = eG.final_CHS_occupancies(tcritical),
        tcritical     = tcritical
      )
)
