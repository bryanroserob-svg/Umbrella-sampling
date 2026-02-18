#!/usr/bin/env python3
"""
jarzynski.py — Estimación Jarzynski ΔG desde SMD (Steered MD)

Calcula el trabajo acumulado W = ∫F·dx usando integración trapezoidal
y estima la energía libre como cota superior (single trajectory).

Uso:
    python3 jarzynski.py <pullf.xvg> <pullx.xvg> <output_dir> <unit>

Argumentos:
    pullf.xvg   - Archivo de fuerza de pulling de GROMACS
    pullx.xvg   - Archivo de posición de pulling de GROMACS
    output_dir   - Directorio de salida para resultados
    unit         - Unidad de energía: 'kCal' o 'kJ'

Salida:
    jarzynski.dat       - Resumen con W total
    jarzynski_work.xvg  - Perfil de trabajo acumulado vs posición
"""

import sys
import os
import numpy as np


def parse_xvg(filepath):
    """Lee un archivo .xvg de GROMACS, ignorando comentarios y metadata."""
    data = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(('#', '@')):
                continue
            vals = line.split()
            if len(vals) >= 2:
                data.append([float(v) for v in vals[:2]])
    return np.array(data) if data else None


def compute_jarzynski(force_data, pos_data):
    """
    Calcula el trabajo acumulado usando integración trapezoidal.

    Parámetros:
        force_data: array Nx2 (tiempo, fuerza en kJ/mol/nm)
        pos_data:   array Nx2 (tiempo, posición en nm)

    Retorna:
        work_cumulative: array de trabajo acumulado
        pos:             array de posiciones correspondientes
        total_work:      trabajo total (kJ/mol)
    """
    n = min(len(force_data), len(pos_data))
    force = force_data[:n, 1]  # kJ/mol/nm
    pos = pos_data[:n, 1]      # nm

    # Trabajo acumulado: W = sum(F_avg * dx)
    dx = np.diff(pos)
    f_avg = (force[:-1] + force[1:]) / 2.0
    work_incremental = f_avg * dx
    work_cumulative = np.cumsum(work_incremental)

    total_work = work_cumulative[-1] if len(work_cumulative) > 0 else 0.0

    return work_cumulative, pos, total_work


def write_results(outdir, work_cumulative, pos, total_work, unit):
    """Escribe los archivos de resultados."""
    os.makedirs(outdir, exist_ok=True)

    factor = 1.0 / 4.184 if 'kCal' in unit else 1.0
    total_work_conv = total_work * factor

    # Resumen
    with open(os.path.join(outdir, 'jarzynski.dat'), 'w') as out:
        out.write("# Jarzynski estimate from SMD\n")
        out.write(f"# W_total = {total_work:.2f} kJ/mol = {total_work/4.184:.2f} kCal/mol\n")
        out.write("# NOTE: Single trajectory — overestimates |ΔG| (upper bound)\n")
        out.write(f"W_kJ {total_work:.4f}\n")
        out.write(f"W_kCal {total_work/4.184:.4f}\n")

    # Perfil de trabajo acumulado
    with open(os.path.join(outdir, 'jarzynski_work.xvg'), 'w') as out:
        out.write("# Cumulative work from SMD\n")
        out.write("# Position(nm)  Work(kJ/mol)\n")
        for i in range(len(work_cumulative)):
            out.write(f"{pos[i+1]:.4f} {work_cumulative[i]:.4f}\n")

    print(f"  ✓ Jarzynski W ≈ {total_work_conv:.1f} {unit}/mol (upper bound)")


def main():
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)

    pullf_file = sys.argv[1]
    pullx_file = sys.argv[2]
    outdir = sys.argv[3]
    unit = sys.argv[4]

    force_data = parse_xvg(pullf_file)
    pos_data = parse_xvg(pullx_file)

    if force_data is None or pos_data is None:
        print("  ⚠ Datos insuficientes para Jarzynski")
        sys.exit(0)

    work_cumulative, pos, total_work = compute_jarzynski(force_data, pos_data)
    write_results(outdir, work_cumulative, pos, total_work, unit)


if __name__ == '__main__':
    main()
