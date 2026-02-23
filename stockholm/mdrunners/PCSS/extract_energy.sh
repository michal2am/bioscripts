#!/bin/bash

# Upewnij się, że moduły są załadowane
# module load gromacs/2026 (lub Twoja wersja)

for sys in sys1 sys2 sys3; do
    echo "Ekstrakcja danych dla $sys..."

    # Przesyłamy numery 15 (Temp), 17 (Pres), 18 (Dens) do gmx energy
    # UWAGA: Numery mogą się różnić zależnie od wersji GROMACS.
    # Możesz sprawdzić poprawne numery wpisując: gmx energy -f sys1/step7_production.edr
    echo "Temperature Pressure Density" | gmx energy -f ${sys}/step7_production.edr -o ${sys}_energy.xvg -b 1000
done