#!/usr/bin/env python3
"""
analyze_convergence.py — Análisis de convergencia por ventana (Block Averaging)

Divide cada trayectoria de umbrella sampling en 5 bloques y calcula:
  - Media y desviación estándar de la posición
  - SEM (Standard Error of the Mean) desde block averaging
  - Drift: diferencia entre el primer y último bloque

Criterios de convergencia:
  - drift < 0.05 nm AND SEM < 0.02 nm → YES
  - De lo contrario → CHECK (necesita más sampling)

Uso:
    python3 analyze_convergence.py <prod_dir> <output_dir> <unit>

Argumentos:
    prod_dir   - Directorio con archivos umbrella*_pullx.xvg
    output_dir - Directorio donde escribir convergence_report.dat
    unit       - Unidad de energía (informativo)
"""

import sys
import os
import glob
import numpy as np


def parse_xvg(filepath):
    """Lee un archivo .xvg de GROMACS, ignorando comentarios y metadata."""
    data = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(('#', '@')):
                continue
            v = line.split()
            if len(v) >= 2:
                data.append([float(x) for x in v[:2]])
    return np.array(data) if data else None


def analyze_window(pullx_file, n_blocks=5):
    """
    Analiza la convergencia de una ventana individual.

    Parámetros:
        pullx_file: path al archivo pullx.xvg de la ventana
        n_blocks:   número de bloques para block averaging

    Retorna:
        tupla (window_name, mean, std, sem, drift, status) o None
    """
    data = parse_xvg(pullx_file)
    if data is None or len(data) < 100:
        return None

    pos = data[:, 1]
    n = len(pos)
    w_name = os.path.basename(pullx_file).replace('_pullx.xvg', '')

    # Block averaging
    block_size = n // n_blocks
    block_means = []
    for b in range(n_blocks):
        block = pos[b * block_size:(b + 1) * block_size]
        block_means.append(np.mean(block))

    overall_mean = np.mean(pos)
    overall_std = np.std(pos)
    block_std = np.std(block_means)
    sem = block_std / np.sqrt(n_blocks)

    # Drift: diferencia entre primer y último bloque
    drift = abs(block_means[-1] - block_means[0])

    converged = "YES" if (drift < 0.05 and sem < 0.02) else "CHECK"
    return (w_name, overall_mean, overall_std, sem, drift, converged)


def write_report(results, output_dir):
    """Escribe el reporte de convergencia."""
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, 'convergence_report.dat')

    bad = 0
    with open(out_file, 'w') as f:
        f.write("# Window  Mean(nm)  Std(nm)  SEM(nm)  Drift(nm)  Status\n")
        for r in results:
            f.write(f"{r[0]:20s} {r[1]:8.4f} {r[2]:8.4f} "
                    f"{r[3]:8.4f} {r[4]:8.4f} {r[5]}\n")
            if r[5] == "CHECK":
                bad += 1

        f.write(f"\n# Total windows: {len(results)}\n")
        f.write(f"# Converged: {len(results) - bad}\n")
        f.write(f"# Need review: {bad}\n")

    print(f"  ✓ Convergencia: {len(results) - bad}/{len(results)} ventanas convergen")
    if bad > 0:
        print(f"  ⚠ {bad} ventanas necesitan más sampling "
              f"(ver convergence_report.dat)")


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    prod_dir = sys.argv[1]
    out_dir = sys.argv[2]
    # unit = sys.argv[3]  # Reserved for future use

    pullx_files = sorted(glob.glob(
        os.path.join(prod_dir, 'umbrella*_pullx.xvg')
    ))

    if not pullx_files:
        print("  ⚠ No se encontraron archivos pullx.xvg")
        sys.exit(0)

    results = []
    for pf in pullx_files:
        result = analyze_window(pf)
        if result is not None:
            results.append(result)

    if results:
        write_report(results, out_dir)
    else:
        print("  ⚠ No hay datos suficientes para análisis de convergencia")


if __name__ == '__main__':
    main()
