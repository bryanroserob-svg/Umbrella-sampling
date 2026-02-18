#!/usr/bin/env python3
"""
plot_umbrella.py — Visualización de resultados de Umbrella Sampling

Genera gráficos de:
  1. Perfil PMF (Potential of Mean Force)
  2. Histogramas de overlap de ventanas
  3. Fuerza de pulling vs tiempo
  4. Distancias COM con ventanas seleccionadas
  5. Convergencia (bootstrap)

Uso:
  python3 plot_umbrella.py <directorio_analisis>
  python3 plot_umbrella.py /ruta/a/US_RUN/.../07_analysis/

Requiere: matplotlib, numpy
"""

import sys
import os
import re
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
except ImportError:
    print("ERROR: matplotlib no está instalado.")
    print("  Instala con: pip install matplotlib")
    sys.exit(1)


# ─── Configuración de estilo ───
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLORS = {
    'pmf': '#2563EB',
    'pmf_fill': '#93C5FD',
    'histogram': None,  # usará colormap
    'force': '#DC2626',
    'distance': '#059669',
    'window_marker': '#F59E0B',
    'bootstrap': '#9333EA',
    'grid': '#E5E7EB',
}


def parse_xvg(filepath):
    """Lee un archivo .xvg de GROMACS y retorna las columnas como arrays numpy."""
    data = []
    metadata = {'title': '', 'xlabel': '', 'ylabel': '', 'legends': []}

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            elif line.startswith('@'):
                # Extraer metadata
                match_title = re.match(r'@\s+title\s+"(.+)"', line)
                match_xaxis = re.match(r'@\s+xaxis\s+label\s+"(.+)"', line)
                match_yaxis = re.match(r'@\s+yaxis\s+label\s+"(.+)"', line)
                match_legend = re.match(r'@\s+s\d+\s+legend\s+"(.+)"', line)

                if match_title:
                    metadata['title'] = match_title.group(1)
                elif match_xaxis:
                    metadata['xlabel'] = match_xaxis.group(1)
                elif match_yaxis:
                    metadata['ylabel'] = match_yaxis.group(1)
                elif match_legend:
                    metadata['legends'].append(match_legend.group(1))
            else:
                try:
                    vals = [float(v) for v in line.split()]
                    if vals:
                        data.append(vals)
                except ValueError:
                    continue

    if not data:
        return None, metadata

    return np.array(data), metadata


def plot_pmf(analysis_dir, output_dir):
    """Genera el gráfico del perfil PMF."""
    pmf_file = os.path.join(analysis_dir, 'pmf.xvg')
    if not os.path.exists(pmf_file):
        print(f"  ⚠ {pmf_file} no encontrado, saltando PMF")
        return

    data, meta = parse_xvg(pmf_file)
    if data is None:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    x = data[:, 0]
    y = data[:, 1]

    # Si hay columna de error (bootstrap)
    if data.shape[1] >= 3:
        yerr = data[:, 2]
        ax.fill_between(x, y - yerr, y + yerr, alpha=0.3, color=COLORS['pmf_fill'],
                        label='Error bootstrap')

    ax.plot(x, y, color=COLORS['pmf'], linewidth=2.5, label='PMF')

    # Anotar ΔG
    min_idx = np.argmin(y)
    max_idx = np.argmax(y)
    delta_g = y[max_idx] - y[min_idx]

    # Encontrar unidad
    unit = 'kCal/mol' if 'kCal' in meta.get('ylabel', '') else 'kJ/mol'
    if 'kCal' not in unit and 'kJ' not in unit:
        unit = 'kCal/mol'

    ax.annotate(f'ΔG ≈ {delta_g:.1f} {unit}',
                xy=(x[max_idx], y[max_idx]),
                xytext=(x[max_idx] - (x[-1]-x[0])*0.15, y[max_idx] * 0.85),
                fontsize=14, fontweight='bold', color=COLORS['pmf'],
                arrowprops=dict(arrowstyle='->', color=COLORS['pmf'], lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=COLORS['pmf'], alpha=0.9))

    ax.set_xlabel(meta.get('xlabel', 'Distancia (nm)'))
    ax.set_ylabel(meta.get('ylabel', f'PMF ({unit})'))
    ax.set_title('Potential of Mean Force (PMF)', fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, color=COLORS['grid'])
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_ylim(bottom=min(y) - abs(delta_g) * 0.1)

    output = os.path.join(output_dir, 'pmf_profile.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ PMF perfil: {output}")

    return delta_g, unit


def plot_histograms(analysis_dir, output_dir):
    """Genera el gráfico de histogramas de overlap."""
    hist_file = os.path.join(analysis_dir, 'histogram.xvg')
    if not os.path.exists(hist_file):
        print(f"  ⚠ {hist_file} no encontrado, saltando histogramas")
        return

    data, meta = parse_xvg(hist_file)
    if data is None:
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    x = data[:, 0]
    n_windows = data.shape[1] - 1

    # Colormap para distinguir ventanas
    cmap = plt.cm.get_cmap('viridis', n_windows)

    for i in range(n_windows):
        color = cmap(i / max(n_windows - 1, 1))
        ax.plot(x, data[:, i + 1], linewidth=1.2, color=color, alpha=0.8)

    ax.set_xlabel(meta.get('xlabel', 'ξ (nm)'))
    ax.set_ylabel(meta.get('ylabel', 'Conteo'))
    ax.set_title(f'Histogramas de Umbrella Sampling ({n_windows} ventanas)', fontweight='bold')
    ax.grid(True, alpha=0.3, color=COLORS['grid'])
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    # Añadir barra de color
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, n_windows - 1))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, label='Ventana #', pad=0.02)
    cbar.set_ticks(np.linspace(0, n_windows - 1, min(n_windows, 10)).astype(int))

    output = os.path.join(output_dir, 'histograms.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ Histogramas: {output}")


def plot_pull_force(analysis_dir, output_dir):
    """Genera el gráfico de fuerza de pulling."""
    force_file = os.path.join(analysis_dir, 'pull_force.xvg')
    if not os.path.exists(force_file):
        print(f"  ⚠ pull_force.xvg no encontrado, saltando")
        return

    data, meta = parse_xvg(force_file)
    if data is None:
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    x = data[:, 0]
    y = data[:, 1]

    ax.plot(x, y, color=COLORS['force'], linewidth=0.5, alpha=0.6, label='Fuerza instantánea')

    # Media móvil para tendencia
    window_size = max(len(x) // 50, 10)
    if len(x) > window_size:
        y_smooth = np.convolve(y, np.ones(window_size) / window_size, mode='valid')
        x_smooth = x[:len(y_smooth)]
        ax.plot(x_smooth, y_smooth, color='#1E3A5F', linewidth=2, label=f'Media móvil (n={window_size})')

    ax.set_xlabel(meta.get('xlabel', 'Tiempo (ps)'))
    ax.set_ylabel(meta.get('ylabel', 'Fuerza (kJ/mol/nm)'))
    ax.set_title('Fuerza de Pulling (SMD)', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, color=COLORS['grid'])
    ax.axhline(y=0, color='gray', linewidth=0.8, linestyle='--')

    output = os.path.join(output_dir, 'pull_force.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ Fuerza de pulling: {output}")


def plot_distances(analysis_dir, output_dir):
    """Genera el gráfico de distancias COM con ventanas marcadas."""
    dist_file = os.path.join(analysis_dir, 'summary_distances.dat')
    windows_file = os.path.join(analysis_dir, 'selected_windows.dat')

    if not os.path.exists(dist_file):
        print(f"  ⚠ summary_distances.dat no encontrado, saltando")
        return

    # Leer distancias
    distances = np.loadtxt(dist_file)
    if distances.ndim != 2 or distances.shape[1] < 2:
        print("  ⚠ Formato inesperado en summary_distances.dat")
        return

    fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(distances[:, 0], distances[:, 1], color=COLORS['distance'],
            linewidth=1.5, label='Distancia COM')

    # Marcar ventanas seleccionadas
    if os.path.exists(windows_file):
        windows = []
        with open(windows_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        windows.append((int(parts[1]), float(parts[2])))
                    except ValueError:
                        continue

        if windows:
            w_frames, w_dists = zip(*windows)
            ax.scatter(w_frames, w_dists, color=COLORS['window_marker'],
                       s=80, zorder=5, edgecolors='black', linewidths=0.8,
                       label=f'Ventanas ({len(windows)})')

    ax.set_xlabel('Frame')
    ax.set_ylabel('Distancia COM (nm)')
    ax.set_title('Distancia COM — Chain A ↔ Chain B', fontweight='bold')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3, color=COLORS['grid'])
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    output = os.path.join(output_dir, 'distances_com.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ Distancias COM: {output}")


def plot_bootstrap(analysis_dir, output_dir):
    """Genera el gráfico de convergencia bootstrap."""
    bs_file = os.path.join(analysis_dir, 'bootstrap_profiles.xvg')
    if not os.path.exists(bs_file):
        print(f"  ⚠ bootstrap_profiles.xvg no encontrado, saltando")
        return

    data, meta = parse_xvg(bs_file)
    if data is None or data.shape[1] < 3:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    x = data[:, 0]
    n_profiles = data.shape[1] - 1

    for i in range(1, min(n_profiles + 1, 51)):
        ax.plot(x, data[:, i], linewidth=0.5, alpha=0.3, color=COLORS['bootstrap'])

    # PMF original superpuesto
    pmf_file = os.path.join(analysis_dir, 'pmf.xvg')
    if os.path.exists(pmf_file):
        pmf_data, _ = parse_xvg(pmf_file)
        if pmf_data is not None:
            ax.plot(pmf_data[:, 0], pmf_data[:, 1], color=COLORS['pmf'],
                    linewidth=2.5, label='PMF', zorder=10)

    ax.set_xlabel(meta.get('xlabel', 'ξ (nm)'))
    ax.set_ylabel(meta.get('ylabel', 'PMF'))
    ax.set_title('Convergencia Bootstrap', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, color=COLORS['grid'])

    output = os.path.join(output_dir, 'bootstrap_convergence.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ Bootstrap: {output}")


def plot_convergence(analysis_dir, output_dir):
    """Genera el gráfico de convergencia por ventana (block averaging)."""
    conv_file = os.path.join(analysis_dir, 'convergence_report.dat')
    if not os.path.exists(conv_file):
        print(f"  ⚠ convergence_report.dat no encontrado, saltando")
        return

    windows, means, stds, sems, drifts, statuses = [], [], [], [], [], []
    with open(conv_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                windows.append(parts[0].replace('umbrella', 'W'))
                means.append(float(parts[1]))
                stds.append(float(parts[2]))
                sems.append(float(parts[3]))
                drifts.append(float(parts[4]))
                statuses.append(parts[5])

    if not windows:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    x = np.arange(len(windows))
    colors_bar = ['#22C55E' if s == 'YES' else '#EF4444' for s in statuses]

    # Panel 1: SEM por ventana
    ax1.bar(x, sems, color=colors_bar, alpha=0.8, edgecolor='white', linewidth=0.5)
    ax1.axhline(y=0.02, color='gray', linewidth=1, linestyle='--', label='Umbral SEM')
    ax1.set_xlabel('Ventana')
    ax1.set_ylabel('SEM (nm)')
    ax1.set_title('Error Estándar por Ventana', fontweight='bold')
    ax1.set_xticks(x[::max(1, len(x)//15)])
    ax1.set_xticklabels([windows[i] for i in range(0, len(windows), max(1, len(x)//15))], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')

    # Panel 2: Drift por ventana
    ax2.bar(x, drifts, color=colors_bar, alpha=0.8, edgecolor='white', linewidth=0.5)
    ax2.axhline(y=0.05, color='gray', linewidth=1, linestyle='--', label='Umbral drift')
    ax2.set_xlabel('Ventana')
    ax2.set_ylabel('Drift (nm)')
    ax2.set_title('Drift por Ventana (Block Avg)', fontweight='bold')
    ax2.set_xticks(x[::max(1, len(x)//15)])
    ax2.set_xticklabels([windows[i] for i in range(0, len(windows), max(1, len(x)//15))], rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    fig.tight_layout()
    output = os.path.join(output_dir, 'convergence.png')
    fig.savefig(output)
    plt.close(fig)

    n_ok = sum(1 for s in statuses if s == 'YES')
    print(f"  ✓ Convergencia: {output} ({n_ok}/{len(statuses)} OK)")


def plot_jarzynski(analysis_dir, output_dir):
    """Genera el gráfico del trabajo acumulado desde SMD (Jarzynski)."""
    work_file = os.path.join(analysis_dir, 'jarzynski_work.xvg')
    if not os.path.exists(work_file):
        print(f"  ⚠ jarzynski_work.xvg no encontrado, saltando")
        return

    data, meta = parse_xvg(work_file)
    if data is None:
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    x = data[:, 0]
    y = data[:, 1]
    y_kcal = y / 4.184

    ax.plot(x, y_kcal, color='#7C3AED', linewidth=2, label='Trabajo acumulado (W)')

    # Comparar con PMF si existe
    pmf_file = os.path.join(analysis_dir, 'pmf.xvg')
    if os.path.exists(pmf_file):
        pmf_data, _ = parse_xvg(pmf_file)
        if pmf_data is not None:
            ax.plot(pmf_data[:, 0], pmf_data[:, 1], color=COLORS['pmf'],
                    linewidth=2, linestyle='--', alpha=0.7, label='PMF (WHAM)')

    ax.set_xlabel('Posición (nm)')
    ax.set_ylabel('W / PMF (kCal/mol)')
    ax.set_title('Jarzynski: Trabajo Acumulado vs PMF', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, color=COLORS['grid'])
    ax.axhline(y=0, color='gray', linewidth=0.8, linestyle=':')

    # Anotar W total
    w_total = y_kcal[-1]
    ax.annotate(f'W = {w_total:.1f} kCal/mol\n(upper bound)',
                xy=(x[-1], y_kcal[-1]), fontsize=11, color='#7C3AED',
                xytext=(-100, 20), textcoords='offset points',
                arrowprops=dict(arrowstyle='->', color='#7C3AED'),
                bbox=dict(boxstyle='round', fc='white', ec='#7C3AED', alpha=0.9))

    output = os.path.join(output_dir, 'jarzynski_work.png')
    fig.savefig(output)
    plt.close(fig)
    print(f"  ✓ Jarzynski: {output}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    analysis_dir = sys.argv[1]
    if not os.path.isdir(analysis_dir):
        print(f"ERROR: {analysis_dir} no es un directorio válido")
        sys.exit(1)

    output_dir = os.path.join(analysis_dir, 'plots')
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n{'='*55}")
    print(f"  UMBRELLA SAMPLING v2 — Visualización")
    print(f"{'='*55}\n")
    print(f"  Directorio: {analysis_dir}")
    print(f"  Salida:     {output_dir}\n")

    # Generar todos los gráficos
    result = plot_pmf(analysis_dir, output_dir)
    plot_histograms(analysis_dir, output_dir)
    plot_pull_force(analysis_dir, output_dir)
    plot_distances(analysis_dir, output_dir)
    plot_bootstrap(analysis_dir, output_dir)
    plot_convergence(analysis_dir, output_dir)    # NUEVO
    plot_jarzynski(analysis_dir, output_dir)      # NUEVO

    print(f"\n{'='*55}")
    if result:
        delta_g, unit = result
        print(f"  ΔG ≈ {delta_g:.2f} {unit}")
    print(f"  Gráficos guardados en: {output_dir}/")
    print(f"{'='*55}\n")


if __name__ == '__main__':
    main()

