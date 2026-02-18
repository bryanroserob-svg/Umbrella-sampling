#!/usr/bin/env python3
"""
validate_overlap.py — Validación cuantitativa de solapamiento entre histogramas

Calcula el integral de solapamiento entre cada par de ventanas adyacentes
de Umbrella Sampling:

    O_ij = ∫ min(P_i(ξ), P_j(ξ)) dξ

donde P_i y P_j son las distribuciones normalizadas de las ventanas i y j.

Criterios:
    - Overlap ≥ 10%  → EXCELENTE (conexión termodinámica robusta)
    - Overlap 5-10%  → OK (aceptable para WHAM)
    - Overlap 1-5%   → POBRE (considerar ventanas intermedias)
    - Overlap < 1%   → CRÍTICO (WHAM puede fallar, agregar ventanas)

Uso:
    python3 validate_overlap.py <histogram.xvg> <output_dir>

Argumentos:
    histogram.xvg - Histogramas generados por gmx wham
    output_dir    - Directorio donde escribir overlap_report.dat

Salida:
    overlap_report.dat  - Reporte detallado por par de ventanas
    Resumen a stdout con calidad global
"""

import sys
import os
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


def compute_overlap(hist_i, hist_j, dx):
    """
    Calcula el integral de solapamiento entre dos histogramas normalizados.

    O_ij = ∫ min(P_i(ξ), P_j(ξ)) dξ

    Parámetros:
        hist_i: array de conteos de la ventana i
        hist_j: array de conteos de la ventana j
        dx:     espaciado en ξ

    Retorna:
        overlap: fracción de solapamiento (0.0 a 1.0)
    """
    # Normalizar a distribuciones de probabilidad
    area_i = np.sum(hist_i) * dx
    area_j = np.sum(hist_j) * dx

    if area_i == 0 or area_j == 0:
        return 0.0

    p_i = hist_i / area_i
    p_j = hist_j / area_j

    # Integral de solapamiento
    overlap = np.sum(np.minimum(p_i, p_j)) * dx

    return overlap


def classify_overlap(overlap_pct):
    """Clasifica la calidad del solapamiento."""
    if overlap_pct >= 10.0:
        return "EXCELENTE"
    elif overlap_pct >= 5.0:
        return "OK"
    elif overlap_pct >= 1.0:
        return "POBRE"
    else:
        return "CRÍTICO"


def analyze_all_overlaps(data):
    """
    Calcula el solapamiento entre todos los pares de ventanas adyacentes.

    Parámetros:
        data: array donde col0=xi, col1..N=histogramas

    Retorna:
        results: lista de tuplas (w_i, w_j, overlap_pct, status)
        global_quality: calidad global del sampling
    """
    xi = data[:, 0]
    hists = data[:, 1:]
    n_windows = hists.shape[1]

    # Espaciado en xi
    dx = xi[1] - xi[0] if len(xi) > 1 else 1.0

    results = []
    for i in range(n_windows - 1):
        overlap = compute_overlap(hists[:, i], hists[:, i + 1], dx)
        overlap_pct = overlap * 100.0
        status = classify_overlap(overlap_pct)
        results.append((i, i + 1, overlap_pct, status))

    # Calidad global
    if not results:
        return results, "N/A"

    min_overlap = min(r[2] for r in results)
    mean_overlap = np.mean([r[2] for r in results])

    if min_overlap >= 5.0:
        global_quality = f"BUENA (min={min_overlap:.1f}%, media={mean_overlap:.1f}%)"
    elif min_overlap >= 1.0:
        global_quality = f"ACEPTABLE (min={min_overlap:.1f}%, media={mean_overlap:.1f}%)"
    else:
        global_quality = f"INSUFICIENTE (min={min_overlap:.1f}%, media={mean_overlap:.1f}%)"

    return results, global_quality


def write_report(results, global_quality, output_dir):
    """Escribe el reporte de solapamiento."""
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, 'overlap_report.dat')

    n_critical = sum(1 for r in results if r[3] == "CRÍTICO")
    n_poor = sum(1 for r in results if r[3] == "POBRE")

    with open(out_file, 'w') as f:
        f.write("# Overlap Validation Report — Umbrella Sampling\n")
        f.write(f"# Global Quality: {global_quality}\n")
        f.write(f"# Total pairs: {len(results)}\n")
        f.write(f"# Critical: {n_critical}, Poor: {n_poor}\n")
        f.write("#\n")
        f.write("# Window_i  Window_j  Overlap(%)  Status\n")
        f.write("# ────────  ────────  ──────────  ──────────\n")
        for w_i, w_j, ovl, status in results:
            f.write(f"  {w_i:>8d}  {w_j:>8d}  {ovl:>10.2f}  {status}\n")

    # Resumen a stdout
    print(f"  ✓ Overlap validation: {len(results)} pares analizados")
    print(f"    Calidad global: {global_quality}")

    if n_critical > 0:
        print(f"    ⚠ {n_critical} pares con overlap CRÍTICO (<1%)")
        print("      → Agregar ventanas intermedias en esas regiones")
        for r in results:
            if r[3] == "CRÍTICO":
                print(f"        Ventanas {r[0]}-{r[1]}: {r[2]:.2f}%")

    if n_poor > 0:
        print(f"    ⚠ {n_poor} pares con overlap POBRE (1-5%)")
        print("      → Considerar ventanas adicionales o más sampling")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    hist_file = sys.argv[1]
    output_dir = sys.argv[2]

    data = parse_histogram(hist_file)
    if data is None:
        print("  ⚠ Histograma vacío, no se puede validar overlap")
        sys.exit(0)

    if data.shape[1] < 3:
        print("  ⚠ Se necesitan al menos 2 ventanas para validar overlap")
        sys.exit(0)

    results, global_quality = analyze_all_overlaps(data)
    write_report(results, global_quality, output_dir)


if __name__ == '__main__':
    main()
