# Umbrella Sampling Pipeline v2 — GROMACS

Automatización completa de umbrella sampling para calcular el **Potential of Mean Force (PMF)** y la energía libre de disociación (**ΔG**).


## 🆕 Mejoras v2

| Feature | Descripción |
|---------|-------------|
| **3 modos** | Multi-cadena, proteína-ligando, permeación de membrana |
| **Carpetas organizadas** | `proteins/`, `ligands/`, `membranes/` con subcarpetas por proyecto |
| **Auto-caja** | Calcula automáticamente centro y dimensiones |
| **Selección adaptativa** | Ventanas por distancia uniforme |
| **Gap detection** | Detecta huecos en histogramas + auto-fill + re-WHAM |
| **Convergencia** | Block averaging por ventana |
| **Jarzynski ΔG** | Cota superior desde SMD |
| **PBC correction** | nojump + center automático |
| **Backup MDP** | `.original` antes de `sed` |
| **Checkpoint granular** | Reanuda ventanas individuales desde `.cpt` |

## Requisitos

- **GROMACS** ≥ 5.x (`gmx` o `gmx_mpi`)
- **Python 3** + `matplotlib` + `numpy`
- **bash** ≥ 4.x, **bc**

## Estructura de carpetas

```
umbrella_sampling/
│
├── umbrella_pipeline.sh          # Script principal (15 etapas)
├── plot_umbrella.py              # 7 tipos de gráficos
├── README.md
│
├── mdp/                          # Archivos MDP (requeridos)
│   ├── ions.mdp
│   ├── minim.mdp
│   ├── npt.mdp
│   ├── md_pull.mdp
│   ├── npt_umbrella.mdp
│   └── md_umbrella.mdp
│
├── proteins/                     # ← Proteínas organizadas por proyecto
│   ├── proteina_tutorial/
│   │   └── 2BEG_model1_capped.pdb
│   ├── caspasa9/
│   │   └── caspasa9_clean.pdb
│   └── receptor_X/
│       └── receptor.gro
│
├── ligands/                      # ← Ligandos parametrizados por proyecto
│   ├── M4_A/
│   │   ├── M4_A.itp             #   Topología (requerido)
│   │   ├── M4_A.gro             #   Coordenadas (recomendado)
│   │   └── M4_A.mol2            #   (opcional, referencia)
│   └── droga_B/
│       ├── droga_B.itp
│       └── droga_B.pdb
│
├── membranes/                    # ← Membranas pre-equilibradas
│   ├── POPC_128/
│   │   └── POPC_128_eq.gro
│   └── DPPC_bilayer/
│       └── DPPC_64.pdb
│
└── US_RUN/                       # (generado automáticamente)
    └── <nombre>_<fecha>/
        ├── 00_setup/ ... 06_umbrella_prod/
        ├── 07_analysis/
        │   ├── pmf.xvg, histogram.xvg
        │   ├── convergence_report.dat
        │   ├── jarzynski.dat
        │   └── plots/           # 7 gráficos PNG
        └── RESUMEN.txt
```

### Menú interactivo de selección

El pipeline muestra las rutas relativas para que identifiques cada archivo:

```
  Proteínas disponibles:
   1) proteins/proteina_tutorial/2BEG_model1_capped.pdb  (76K)
   2) proteins/caspasa9/caspasa9_clean.pdb  (245K)
  Seleccione proteína [1-2]: █

  Ligandos parametrizados:
   1) ligands/M4_A/M4_A.itp  + .gro
  Seleccione ligando [1-1] (0=omitir): █
```

## Cómo preparar cada carpeta

### proteins/
Crea una **subcarpeta por proteína** con el PDB o GRO dentro:
```bash
mkdir -p proteins/mi_proteina/
cp mi_proteina.pdb proteins/mi_proteina/
```

### ligands/ (debe estar parametrizado)
Los ligandos **necesitan topología GROMACS** (`.itp`) antes de usarlos. Crea una subcarpeta por ligando:

```bash
mkdir -p ligands/mi_ligando/
# Copia el .itp + coordenadas:
cp mi_ligando.itp ligands/mi_ligando/
cp mi_ligando.gro ligands/mi_ligando/   # o .pdb
```

**Cómo parametrizar:**

| Herramienta | Force Field | Uso |
|-------------|------------|-----|
| **CGenFF** | CHARMM | [cgenff.silcsbio.org](https://cgenff.silcsbio.org/) — sube `.mol2` |
| **ACPYPE** | AMBER/GAFF | `acpype -i ligando.mol2 -c bcc` |
| **SwissParam** | CHARMM | [swissparam.ch](http://swissparam.ch/) |
| **ATB** | GROMOS | [atb.uq.edu.au](https://atb.uq.edu.au/) |

### membranes/
Necesitas una membrana **pre-equilibrada**. Subcarpeta por tipo:
```bash
mkdir -p membranes/POPC_128/
# Genera con CHARMM-GUI, insane.py, o PACKMOL
```

## Uso

### Ejecución

```bash
chmod +x umbrella_pipeline.sh
./umbrella_pipeline.sh
```

### Reanudar (checkpoint)

```bash
./umbrella_pipeline.sh --resume US_RUN/<nombre>/
# Reanuda etapas completas Y ventanas individuales desde .cpt
```

### Graficar

```bash
python3 plot_umbrella.py US_RUN/<nombre>/07_analysis/
```

## Modos del sistema

| Modo | Grupo móvil | Grupo referencia | Uso |
|------|------------|-----------------|-----|
| **multichain** | Chain_A (residuos) | Chain_B (residuos) | Separar cadenas |
| **protlig** | Ligando | Proteína | Unbinding |
| **membrane** | Permeant | Membrana | Permeación |

## Etapas del pipeline (15)

| # | Etapa | Descripción |
|---|-------|-------------|
| 1 | Topología | `gmx pdb2gmx` |
| 2 | Restraints | POSRES_B (solo multichain) |
| 3 | Caja | Auto-detect o manual |
| 4 | Solvatar + iones | `gmx solvate` + `gmx genion` |
| 5 | Minimización | Steepest descent |
| 6 | NPT eq. | Equilibración NPT |
| 7 | Index groups | Según modo |
| 8 | Pulling SMD | Steered MD |
| 8b | PBC correction | nojump + center |
| 9 | Frames + distancias | Extraer + COM dist |
| 9b | Jarzynski | W = ∫F·dx |
| 10 | Ventanas adaptativas | Greedy por distancia |
| 11-12 | Umbrella NPT + prod | Checkpoint granular |
| 13 | WHAM | PMF + bootstrap |
| 14 | Convergencia | Block averaging |
| 15 | Gap detection | Auto-fill + re-WHAM |

## Gráficos generados (7)

| Archivo | Contenido |
|---------|-----------|
| `pmf_profile.png` | Perfil PMF con ΔG |
| `histograms.png` | Overlap de histogramas |
| `pull_force.png` | Fuerza SMD |
| `distances_com.png` | Distancias COM + ventanas |
| `bootstrap_convergence.png` | Perfiles bootstrap |
| `convergence.png` | SEM + drift por ventana |
| `jarzynski_work.png` | Trabajo acumulado vs PMF |

## Troubleshooting

| Error | Solución |
|-------|----------|
| No encuentra PDBs | Crea `proteins/<nombre>/tu_proteina.pdb` |
| Ligando no encontrado | Parametrizar → copiar `.itp` + `.gro` a `ligands/<nombre>/` |
| Histogramas sin overlap | Reducir `WINDOW_SPACING` o aumentar `PULL_K` |
| WHAM no converge | Más tiempo de producción |
| Ventanas `CHECK` | Extender esas ventanas |
| Ventana interrumpida | `--resume` detecta `.cpt` y continúa |
