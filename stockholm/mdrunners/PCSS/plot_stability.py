import matplotlib.pyplot as plt
import numpy as np
import os


def read_xvg(filename):
    """Funkcja do czytania plików XVG z pominięciem nagłówków GROMACS."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                data.append([float(x) for x in line.split()])
    return np.array(data)


systems = ['sys1', 'sys2', 'sys3']
metrics = ['Temperature [K]', 'Pressure [bar]', 'Density [kg/m^3]']

fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
plt.subplots_adjust(hspace=0.3)

for sys_name in systems:
    fname = f"{sys_name}_energy.xvg"
    if not os.path.exists(fname):
        print(f"Ominięto {sys_name} - brak pliku {fname}")
        continue

    data = read_xvg(fname)
    time = data[:, 0] / 1000  # Zamiana ps na ns

    for i in range(3):
        # i+1 bo kolumna 0 to czas
        axs[i].plot(time, data[:, i + 1], label=f"{sys_name}", alpha=0.7)
        axs[i].set_ylabel(metrics[i])
        axs[i].grid(True, linestyle='--', alpha=0.5)

axs[2].set_xlabel("Time [ns]")
axs[0].legend(loc='upper right', ncol=3)
axs[0].set_title("Analiza stabilności systemów (GABA + Membrana)")

plt.savefig('stability_check.png', dpi=300, bbox_inches='tight')
print("Sukces! Wykres został zapisany jako: stability_check.png")