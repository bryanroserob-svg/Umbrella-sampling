#!/usr/bin/env python3
"""
select_windows.py — Selección adaptativa de ventanas para Umbrella Sampling

Selecciona frames cuyas distancias COM estén uniformemente espaciadas
usando un algoritmo greedy con tolerancia del 90% del espaciado objetivo.

Uso:
    python3 select_windows.py <spacing_nm> <summary_distances.dat> <output_file>

Argumentos:
    spacing_nm            - Espaciado deseado entre ventanas (nm)
    summary_distances.dat - Archivo con columnas: frame distancia
    output_file           - Archivo de salida: columnas window frame distance

Salida:
    Imprime WINDOWS=N y RANGE=min-max a stdout para captura por bash.
"""

import sys
import numpy as np


def select_adaptive_windows(distances_file, spacing):
    """
    Selecciona ventanas con espaciado uniforme usando algoritmo greedy.

    Parámetros:
        distances_file: path al archivo summary_distances.dat
        spacing:        espaciado objetivo en nm

    Retorna:
        list of tuples (frame, distance)
    """
    data = np.loadtxt(distances_file)
    frames = data[:, 0].astype(int)
    distances = data[:, 1]

    # Ordenar por distancia
    idx_sorted = np.argsort(distances)
    frames_sorted = frames[idx_sorted]
    dist_sorted = distances[idx_sorted]

    # Seleccionar ventanas: greedy, cada ~spacing nm
    selected = [(frames_sorted[0], dist_sorted[0])]
    last_d = dist_sorted[0]

    for i in range(1, len(dist_sorted)):
        if abs(dist_sorted[i] - last_d) >= spacing * 0.9:  # 90% del spacing
            selected.append((frames_sorted[i], dist_sorted[i]))
            last_d = dist_sorted[i]

    return selected


def write_output(selected, output_file):
    """Escribe el archivo de ventanas seleccionadas."""
    with open(output_file, 'w') as f:
        f.write("# window frame distance\n")
        for w, (fr, d) in enumerate(selected):
            f.write(f"{w} {fr} {d:.4f}\n")

    # Estos prints son capturados por bash
    print(f"WINDOWS={len(selected)}")
    print(f"RANGE={selected[0][1]:.3f}-{selected[-1][1]:.3f}")


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    spacing = float(sys.argv[1])
    dist_file = sys.argv[2]
    out_file = sys.argv[3]

    selected = select_adaptive_windows(dist_file, spacing)
    write_output(selected, out_file)


if __name__ == '__main__':
    main()
