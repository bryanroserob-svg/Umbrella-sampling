#!/usr/bin/env python3
"""
detect_gaps.py — Detección de gaps en histogramas de Umbrella Sampling

Analiza el archivo histogram.xvg generado por gmx wham para detectar
regiones de ξ donde ninguna ventana tiene muestreo significativo.
Para cada gap, busca el frame más cercano en summary_distances.dat
para sugerir ventanas de relleno (gap-filling).

Uso:
    python3 detect_gaps.py <histogram.xvg> <spacing> <summary_distances.dat> <gap_output>

Argumentos:
    histogram.xvg         - Histogramas de WHAM
    spacing               - Espaciado de ventanas (nm)
    summary_distances.dat - Distancias COM por frame
    gap_output            - Archivo de salida con ventanas sugeridas
"""

import sys
import numpy as np


def parse_histogram(hist_file):
    """Lee el archivo de histogramas de WHAM."""
    data = []
    with open(hist_file) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            vals = line.split()
            if vals:
                data.append([float(v) for v in vals])
    return np.array(data) if data else None


def find_gaps(data, spacing):
    """
    Detecta regiones sin cobertura en los histogramas.

    Parámetros:
        data:    array donde col0=xi, col1..N=histogramas
        spacing: espaciado de ventanas para agrupar gaps

    Retorna:
        lista de centros de gap (nm)
    """
    xi = data[:, 0]
    hists = data[:, 1:]

    # Umbral: 1% del máximo promedio de cada ventana
    max_per_window = np.max(hists, axis=0)
    threshold = np.mean(max_per_window) * 0.01

    # Puntos xi donde la suma total es menor al umbral
    total_count = np.sum(hists, axis=1)
    gap_regions = xi[total_count < threshold]

    if len(gap_regions) == 0:
        return []

    # Agrupar regiones contiguas y calcular centros
    gap_centers = []
    current_gap = [gap_regions[0]]
    for i in range(1, len(gap_regions)):
        if gap_regions[i] - gap_regions[i - 1] < spacing:
            current_gap.append(gap_regions[i])
        else:
            gap_centers.append(np.mean(current_gap))
            current_gap = [gap_regions[i]]
    gap_centers.append(np.mean(current_gap))

    return gap_centers


def find_nearest_frames(gap_centers, dist_file):
    """
    Para cada centro de gap, encuentra el frame más cercano
    en summary_distances.dat.

    Retorna:
        list of tuples (frame, actual_distance, gap_center)
    """
    dist_data = np.loadtxt(dist_file)
    frames_all = dist_data[:, 0].astype(int)
    dists_all = dist_data[:, 1]

    fill_windows = []
    for gc in gap_centers:
        idx = np.argmin(np.abs(dists_all - gc))
        fill_windows.append((frames_all[idx], dists_all[idx], gc))

    return fill_windows


def write_results(fill_windows, gap_centers, gap_file):
    """Escribe las ventanas sugeridas para rellenar gaps."""
    if not gap_centers:
        print("  ✓ No se detectaron gaps en los histogramas")
        with open(gap_file, 'w') as f:
            f.write("# No gaps detected\n")
        return

    with open(gap_file, 'w') as f:
        f.write("# Gap-filling windows\n")
        f.write("# frame  actual_dist  gap_center\n")
        for fw, ad, gc in fill_windows:
            f.write(f"{fw} {ad:.4f} {gc:.4f}\n")

    print(f"  ⚠ Se detectaron {len(gap_centers)} gaps en los histogramas")
    for gc in gap_centers:
        print(f"    Gap en ξ ≈ {gc:.3f} nm")
    print(f"  → Ventanas sugeridas guardadas en {gap_file}")


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    hist_file = sys.argv[1]
    spacing = float(sys.argv[2])
    dist_file = sys.argv[3]
    gap_file = sys.argv[4]

    data = parse_histogram(hist_file)
    if data is None:
        print("  ⚠ Histograma vacío")
        sys.exit(0)

    gap_centers = find_gaps(data, spacing)

    if gap_centers:
        fill_windows = find_nearest_frames(gap_centers, dist_file)
    else:
        fill_windows = []

    write_results(fill_windows, gap_centers, gap_file)


if __name__ == '__main__':
    main()
