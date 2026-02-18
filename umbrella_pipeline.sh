#!/bin/bash
set -euo pipefail

###############################################################################
#  UMBRELLA SAMPLING PIPELINE v2 — GROMACS
#
#  Automatización completa del flujo de trabajo de umbrella sampling.
#  Soporta 3 modos: multi-chain, proteína-ligando, permeación membrana.
#
#  Mejoras v2:
#    - Modo proteína-ligando y permeación de membrana
#    - Auto-detección de dimensiones de caja
#    - Selección adaptativa de ventanas (por distancia)
#    - Detección de gaps en histogramas + auto-fill
#    - Análisis de convergencia por ventana (block averaging)
#    - Estimación Jarzynski desde SMD
#    - Corrección automática de PBC
#
#  Requiere: GROMACS >= 5.x, bash >= 4.x
###############################################################################

#==========================================
# CONFIGURACIÓN
#==========================================
if command -v gmx &> /dev/null; then
    readonly GMX=gmx; readonly USE_MPI=false
elif command -v gmx_mpi &> /dev/null; then
    readonly GMX=gmx_mpi; readonly USE_MPI=true
else
    echo "ERROR: No se encontró GROMACS (gmx o gmx_mpi)" >&2; exit 1
fi

readonly NT=${UMBRELLA_NT:-$(nproc --all 2>/dev/null || echo 4)}
if [ "$USE_MPI" = true ]; then
    readonly MDRUN="$GMX mdrun"
else
    readonly MDRUN="$GMX mdrun -nt $NT"
fi

readonly BASE_MDP=mdp
readonly INITIAL_DIR="$(pwd)"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
WORKDIR="US_RUN"

# Modo del sistema (se define en get_user_input)
SYSTEM_MODE=""  # multichain | protlig | membrane

#==========================================
# COLORES
#==========================================
readonly RED='\033[0;31m' GREEN='\033[0;32m' YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m' CYAN='\033[0;36m' MAGENTA='\033[0;35m' NC='\033[0m'
CURRENT_STAGE="inicialización"

#==========================================
# FUNCIONES AUXILIARES
#==========================================
log_step() {
    CURRENT_STAGE="$1"
    echo -e "\n${BLUE}=========================================${NC}"
    echo -e "${BLUE}>>> $1${NC}"
    echo -e "${BLUE}=========================================${NC}\n"
}
log_success() { echo -e "${GREEN}✓${NC} $1"; }
log_error()   { echo -e "${RED}❌ ERROR:${NC} $1" >&2; }
log_warning() { echo -e "${YELLOW}⚠${NC} $1"; }
log_info()    { echo -e "${CYAN}ℹ${NC} $1"; }

create_dir() { [ ! -d "$1" ] && mkdir -p "$1" && log_success "Carpeta: $1"; }
run_gmx() { "$GMX" "$@"; }

mark_stage_done() { echo "$1" >> "$RUNDIR/.completed_stages"; }
is_stage_done() { [ -f "$RUNDIR/.completed_stages" ] && grep -qxF "$1" "$RUNDIR/.completed_stages"; }

# Backup de MDPs antes de modificar con sed
backup_mdp() {
    local mdp_file="$1"
    if [ -f "$mdp_file" ] && [ ! -f "${mdp_file}.original" ]; then
        cp "$mdp_file" "${mdp_file}.original"
    fi
}

# Carpetas de entrada (se crean si no existen)
readonly INPUT_PROTEINS_DIR="${INITIAL_DIR}/proteins"
readonly INPUT_LIGANDS_DIR="${INITIAL_DIR}/ligands"
readonly INPUT_MEMBRANES_DIR="${INITIAL_DIR}/membranes"

cleanup_on_error() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo ""; log_error "Interrumpido durante: ${CURRENT_STAGE}"
        log_info "Reanudar con: $0 --resume ${RUNDIR:-}"
    fi
}
trap cleanup_on_error EXIT

#==========================================
# AYUDA
#==========================================
show_help() {
    cat <<EOF
${BLUE}  UMBRELLA SAMPLING PIPELINE v3 — GROMACS${NC}

${CYAN}Uso:${NC}
  $0                    Ejecución interactiva
  $0 --resume <DIR>     Reanudar
  $0 --cleanup <DIR>    Limpiar archivos temporales
  $0 --help             Ayuda

${CYAN}Modos:${NC}
  1) Multi-cadena   — Separar cadenas proteicas (tutorial Lemkul)
  2) Proteína-Ligando — Extraer ligando del sitio activo
  3) Permeación      — Molécula a través de membrana lipídica

${CYAN}Features v3:${NC}
  • Scripts Python modulares (scripts/)
  • Selección adaptativa de ventanas por distancia
  • Detección de gaps en histogramas + auto-fill
  • Validación cuantitativa de overlap entre ventanas
  • Convergencia por ventana (block averaging)
  • Estimación Jarzynski ΔG desde SMD
  • Corrección automática de PBC
  • Auto-detección de dimensiones de caja
  • Limpieza automática de temporales (--cleanup)
EOF
    exit 0
}

#==========================================
# INPUT DEL USUARIO
#==========================================
get_user_input() {
    echo -e "${MAGENTA}═══════════════════════════════════════════${NC}"
    echo -e "${MAGENTA}  UMBRELLA SAMPLING v2 — CONFIGURACIÓN${NC}"
    echo -e "${MAGENTA}═══════════════════════════════════════════${NC}\n"

    # ── Modo del sistema ──
    echo -e "${CYAN}═══ MODO DEL SISTEMA ═══${NC}\n"
    echo "   1) Multi-cadena    — Separar una cadena de otra (tutorial Lemkul)"
    echo "   2) Proteína-Ligando — Extraer un ligando del sitio activo"
    echo "   3) Permeación      — Molécula a través de membrana lipídica"
    echo -e "\n  ${YELLOW}Seleccione modo [1-3] (default: 1):${NC}"
    read -r MODE_CHOICE
    case $MODE_CHOICE in
        2) SYSTEM_MODE="protlig" ;;
        3) SYSTEM_MODE="membrane" ;;
        *) SYSTEM_MODE="multichain" ;;
    esac
    log_success "Modo: $SYSTEM_MODE"

    # ── Selección de proteína (desde carpeta proteins/ con subcarpetas) ──
    echo -e "\n${CYAN}═══ ESTRUCTURA PROTEICA ═══${NC}\n"
    local pdb_files=()

    # Buscar recursivamente en proteins/
    if [ -d "$INPUT_PROTEINS_DIR" ]; then
        while IFS= read -r -d '' f; do
            pdb_files+=("$f")
        done < <(find "$INPUT_PROTEINS_DIR" \( -name "*.pdb" -o -name "*.gro" \) -type f -print0 | sort -z)
        if [ ${#pdb_files[@]} -gt 0 ]; then
            log_info "Archivos encontrados en proteins/"
        fi
    fi

    # Fallback: buscar en directorio actual
    if [ ${#pdb_files[@]} -eq 0 ]; then
        while IFS= read -r -d '' f; do
            pdb_files+=("$f")
        done < <(find "$INITIAL_DIR" -maxdepth 1 -name "*.pdb" -type f -print0 | sort -z)
    fi

    if [ ${#pdb_files[@]} -eq 0 ]; then
        log_error "No se encontraron archivos .pdb ni en proteins/ ni en el directorio actual"
        echo -e "  ${CYAN}Crea la carpeta proteins/<nombre_proyecto>/:${NC}"
        echo "    mkdir -p proteins/mi_proteina/"
        echo "    cp tu_proteina.pdb proteins/mi_proteina/"
        exit 1
    fi

    echo "  Proteínas disponibles:"
    for i in "${!pdb_files[@]}"; do
        local fsize rel_path
        fsize=$(du -h "${pdb_files[$i]}" | cut -f1)
        rel_path=$(realpath --relative-to="$INITIAL_DIR" "${pdb_files[$i]}")
        echo "   $((i+1))) $rel_path  ($fsize)"
    done
    echo -e "\n  ${YELLOW}Seleccione proteína [1-${#pdb_files[@]}]:${NC}"
    read -r PDB_CHOICE
    if [[ "$PDB_CHOICE" =~ ^[0-9]+$ ]] && [ "$PDB_CHOICE" -ge 1 ] && [ "$PDB_CHOICE" -le ${#pdb_files[@]} ]; then
        INPUT_PDB="${pdb_files[$((PDB_CHOICE-1))]}"
    else
        log_error "Selección inválida"; exit 1
    fi
    log_success "Proteína: $(basename "$INPUT_PDB") [$(dirname "$(realpath --relative-to="$INITIAL_DIR" "$INPUT_PDB")")]"

    # ── Selección de ligando (solo modo protlig/membrane) ──
    INPUT_LIGAND_ITP=""
    INPUT_LIGAND_GRO=""
    if [ "$SYSTEM_MODE" = "protlig" ] || [ "$SYSTEM_MODE" = "membrane" ]; then
        echo -e "\n${CYAN}═══ LIGANDO / MOLÉCULA ═══${NC}\n"
        if [ -d "$INPUT_LIGANDS_DIR" ]; then
            local lig_files=()
            while IFS= read -r -d '' f; do
                lig_files+=("$f")
            done < <(find "$INPUT_LIGANDS_DIR" -name "*.itp" -type f -print0 | sort -z)

            if [ ${#lig_files[@]} -gt 0 ]; then
                echo -e "  ${GREEN}Ligandos parametrizados:${NC}"
                for i in "${!lig_files[@]}"; do
                    local lig_dir lig_base has_coords
                    lig_base=$(basename "${lig_files[$i]}" .itp)
                    lig_dir=$(dirname "${lig_files[$i]}")
                    has_coords=""
                    [ -f "$lig_dir/${lig_base}.gro" ] && has_coords="+ .gro"
                    [ -f "$lig_dir/${lig_base}.pdb" ] && has_coords="+ .pdb"
                    local rel_path
                    rel_path=$(realpath --relative-to="$INITIAL_DIR" "${lig_files[$i]}")
                    echo "   $((i+1))) $rel_path  $has_coords"
                done
                echo -e "\n  ${YELLOW}Seleccione ligando [1-${#lig_files[@]}] (0=omitir):${NC}"
                read -r LIG_CHOICE
                if [[ "$LIG_CHOICE" =~ ^[1-9][0-9]*$ ]] && [ "$LIG_CHOICE" -le ${#lig_files[@]} ]; then
                    INPUT_LIGAND_ITP="${lig_files[$((LIG_CHOICE-1))]}"
                    local lig_base lig_dir
                    lig_base=$(basename "$INPUT_LIGAND_ITP" .itp)
                    lig_dir=$(dirname "$INPUT_LIGAND_ITP")
                    [ -f "$lig_dir/${lig_base}.gro" ] && INPUT_LIGAND_GRO="$lig_dir/${lig_base}.gro"
                    [ -f "$lig_dir/${lig_base}.pdb" ] && INPUT_LIGAND_GRO="$lig_dir/${lig_base}.pdb"
                    log_success "Ligando: $lig_base (.itp${INPUT_LIGAND_GRO:+ + coord})"
                fi
            else
                log_warning "No se encontraron .itp en ligands/"
                echo -e "  ${CYAN}Para usar ligandos, parametrízalos primero:${NC}"
                echo "    1. Genera topología con CGenFF, ACPYPE, o ATB"
                echo "    2. Crea subcarpeta: mkdir -p ligands/mi_ligando/"
                echo "    3. Copia el .itp y .gro/.pdb a la subcarpeta"
                echo "    4. Re-ejecuta el pipeline"
            fi
        else
            log_warning "Carpeta ligands/ no encontrada"
            echo -e "  ${CYAN}Crea ligands/ con subcarpetas para cada ligando:${NC}"
            echo "    mkdir -p ligands/mi_ligando/"
        fi
    fi

    # ── Selección de membrana (solo modo membrane) ──
    INPUT_MEMBRANE=""
    if [ "$SYSTEM_MODE" = "membrane" ]; then
        echo -e "\n${CYAN}═══ MEMBRANA ═══${NC}\n"
        if [ -d "$INPUT_MEMBRANES_DIR" ]; then
            local mem_files=()
            while IFS= read -r -d '' f; do
                mem_files+=("$f")
            done < <(find "$INPUT_MEMBRANES_DIR" \( -name "*.pdb" -o -name "*.gro" \) -type f -print0 | sort -z)

            if [ ${#mem_files[@]} -gt 0 ]; then
                echo "  Membranas disponibles:"
                for i in "${!mem_files[@]}"; do
                    local rel_path
                    rel_path=$(realpath --relative-to="$INITIAL_DIR" "${mem_files[$i]}")
                    echo "   $((i+1))) $rel_path"
                done
                echo -e "\n  ${YELLOW}Seleccione membrana [1-${#mem_files[@]}] (0=omitir):${NC}"
                read -r MEM_CHOICE
                if [[ "$MEM_CHOICE" =~ ^[1-9][0-9]*$ ]] && [ "$MEM_CHOICE" -le ${#mem_files[@]} ]; then
                    INPUT_MEMBRANE="${mem_files[$((MEM_CHOICE-1))]}"
                    log_success "Membrana: $(basename "$INPUT_MEMBRANE")"
                fi
            fi
        else
            log_warning "Carpeta membranes/ no encontrada"
            echo -e "  ${CYAN}Crea membranes/ con subcarpetas:${NC}"
            echo "    mkdir -p membranes/POPC_128/"
            echo "    # Genera con CHARMM-GUI o insane.py"
        fi
    fi

    # ── Force field ──
    echo -e "\n${CYAN}═══ FORCE FIELD ═══${NC}\n"
    echo -e "  ${YELLOW}Número de force field para pdb2gmx (enter=interactivo):${NC}"
    read -r FF_NUMBER; FF_NUMBER=${FF_NUMBER:-""}

    # ── Modelo de agua ──
    echo -e "\n${CYAN}═══ MODELO DE AGUA ═══${NC}\n"
    echo "   1) spc216  2) tip3p  3) tip4p"
    echo -e "  ${YELLOW}[1-3] (default: 1):${NC}"
    read -r WATER_CHOICE
    case $WATER_CHOICE in
        2) WATER_MODEL="tip3p"; WATER_FILE="spc216.gro" ;;
        3) WATER_MODEL="tip4p"; WATER_FILE="tip4p.gro" ;;
        *) WATER_MODEL="spc"; WATER_FILE="spc216.gro" ;;
    esac

    # ── Grupos de pulling según modo ──
    echo -e "\n${CYAN}═══ GRUPOS DE PULLING ═══${NC}\n"
    case $SYSTEM_MODE in
        multichain)
            echo -e "  ${YELLOW}Residuos grupo MOVIBLE (Chain_A), ej: 1-27:${NC}"
            read -r CHAIN_A_RESIDUES; CHAIN_A_RESIDUES=${CHAIN_A_RESIDUES:-"1-27"}
            echo -e "  ${YELLOW}Residuos grupo REFERENCIA (Chain_B), ej: 28-54:${NC}"
            read -r CHAIN_B_RESIDUES; CHAIN_B_RESIDUES=${CHAIN_B_RESIDUES:-"28-54"}
            echo -e "  ${YELLOW}Nombre .itp cadena referencia (ej: Protein_chain_B):${NC}"
            read -r REF_CHAIN_ITP; REF_CHAIN_ITP=${REF_CHAIN_ITP:-"Protein_chain_B"}
            PULL_GROUP1_NAME="Chain_A"; PULL_GROUP2_NAME="Chain_B"
            ;;
        protlig)
            echo -e "  ${CYAN}Grupo móvil = Ligando, Grupo referencia = Proteína${NC}"
            echo -e "  ${YELLOW}Nombre del residuo del ligando [default: LIG]:${NC}"
            read -r LIG_RESNAME; LIG_RESNAME=${LIG_RESNAME:-"LIG"}
            CHAIN_A_RESIDUES="$LIG_RESNAME"; CHAIN_B_RESIDUES="Protein"
            REF_CHAIN_ITP=""
            PULL_GROUP1_NAME="LIG"; PULL_GROUP2_NAME="Protein"
            ;;
        membrane)
            echo -e "  ${CYAN}Grupo móvil = molécula a permear, Grupo ref = membrana COM${NC}"
            echo -e "  ${YELLOW}Residuos/nombre de la molécula a permear [ej: LIG o 1-10]:${NC}"
            read -r CHAIN_A_RESIDUES; CHAIN_A_RESIDUES=${CHAIN_A_RESIDUES:-"LIG"}
            echo -e "  ${YELLOW}Nombre del lípido de la membrana [default: POPC]:${NC}"
            read -r MEMBRANE_LIPID; MEMBRANE_LIPID=${MEMBRANE_LIPID:-"POPC"}
            CHAIN_B_RESIDUES="$MEMBRANE_LIPID"
            REF_CHAIN_ITP=""
            PULL_GROUP1_NAME="Permeant"; PULL_GROUP2_NAME="Membrane"
            ;;
    esac
    log_success "Grupo móvil: $CHAIN_A_RESIDUES | Referencia: $CHAIN_B_RESIDUES"

    # ── Caja de simulación ──
    echo -e "\n${CYAN}═══ CAJA DE SIMULACIÓN ═══${NC}\n"
    echo "   1) Auto-detectar dimensiones (recomendado)"
    echo "   2) Especificar manualmente"
    echo -e "  ${YELLOW}[1-2] (default: 1):${NC}"
    read -r BOX_MODE_CHOICE
    BOX_MODE=${BOX_MODE_CHOICE:-1}

    if [ "$BOX_MODE" = "2" ]; then
        echo -e "  ${YELLOW}Centro (x y z nm), ej: 3.280 2.181 2.4775:${NC}"
        read -r BOX_CENTER; BOX_CENTER=${BOX_CENTER:-"3.280 2.181 2.4775"}
        echo -e "  ${YELLOW}Dimensiones (x y z nm), ej: 6.560 4.362 12.000:${NC}"
        read -r BOX_DIMS; BOX_DIMS=${BOX_DIMS:-"6.560 4.362 12.000"}
    else
        BOX_CENTER="auto"; BOX_DIMS="auto"
        echo -e "  ${YELLOW}Distancia de pulling deseada (nm) [default: 5.0]:${NC}"
        read -r PULL_DISTANCE; PULL_DISTANCE=${PULL_DISTANCE:-5.0}
        echo -e "  ${YELLOW}Margen extra en ejes no-pulling (nm) [default: 1.5]:${NC}"
        read -r BOX_MARGIN; BOX_MARGIN=${BOX_MARGIN:-1.5}
    fi

    # ── Dirección de pulling ──
    echo -e "\n${CYAN}═══ DIRECCIÓN DE PULLING ═══${NC}\n"
    echo "   1) Solo Z (N N Y)  2) Solo X (Y N N)  3) Solo Y (N Y N)  4) 3D (Y Y Y)"
    echo -e "  ${YELLOW}[1-4] (default: 1):${NC}"
    read -r PULL_DIR_CHOICE
    case $PULL_DIR_CHOICE in
        2) PULL_DIM="Y N N"; PULL_AXIS="x" ;;
        3) PULL_DIM="N Y N"; PULL_AXIS="y" ;;
        4) PULL_DIM="Y Y Y"; PULL_AXIS="xyz" ;;
        *) PULL_DIM="N N Y"; PULL_AXIS="z" ;;
    esac

    # ── Iones ──
    echo -e "\n${CYAN}═══ IONES ═══${NC}"
    echo -e "  ${YELLOW}Concentración NaCl (mol/L) [default: 0.1]:${NC}"
    read -r ION_CONC; ION_CONC=${ION_CONC:-0.1}

    # ── Parámetros de pulling ──
    echo -e "\n${CYAN}═══ PULLING (SMD) ═══${NC}"
    echo -e "  ${YELLOW}Velocidad (nm/ps) [0.01]:${NC}"; read -r PULL_RATE; PULL_RATE=${PULL_RATE:-0.01}
    echo -e "  ${YELLOW}Constante k (kJ/mol/nm²) [1000]:${NC}"; read -r PULL_K; PULL_K=${PULL_K:-1000}
    echo -e "  ${YELLOW}Tiempo pulling (ps) [500]:${NC}"; read -r PULL_TIME; PULL_TIME=${PULL_TIME:-500}
    PULL_NSTEPS=$(echo "$PULL_TIME / 0.002" | bc)

    # ── Parámetros umbrella ──
    echo -e "\n${CYAN}═══ UMBRELLA SAMPLING ═══${NC}"
    echo -e "  ${YELLOW}Espaciado entre ventanas (nm) [0.2]:${NC}"; read -r WINDOW_SPACING; WINDOW_SPACING=${WINDOW_SPACING:-0.2}
    echo -e "  ${YELLOW}Tiempo producción (ns) [10]:${NC}"; read -r UMBRELLA_NS; UMBRELLA_NS=${UMBRELLA_NS:-10}
    UMBRELLA_NSTEPS=$(echo "$UMBRELLA_NS * 500000" | bc)
    echo -e "  ${YELLOW}Unidades PMF (kCal/kJ) [kCal]:${NC}"; read -r PMF_UNIT; PMF_UNIT=${PMF_UNIT:-kCal}

    PARALLEL_WINDOWS=${UMBRELLA_PARALLEL:-1}
    echo -e "  ${YELLOW}Ventanas en paralelo [${PARALLEL_WINDOWS}]:${NC}"; read -r p; PARALLEL_WINDOWS=${p:-$PARALLEL_WINDOWS}

    # Guardar el frame increment como fallback
    FRAME_INCREMENT=20

    log_success "Configuración completada"
}

#==========================================
# VALIDACIÓN MDP
#==========================================
validate_mdp_files() {
    log_step "Verificando archivos MDP"
    [ ! -d "$BASE_MDP" ] && { log_error "Carpeta '$BASE_MDP/' no encontrada"; exit 1; }
    for mdp in ions.mdp minim.mdp npt.mdp md_pull.mdp npt_umbrella.mdp md_umbrella.mdp; do
        [ ! -f "$BASE_MDP/$mdp" ] && { log_error "$BASE_MDP/$mdp no encontrado"; exit 1; }
        log_success "$mdp ✓"
    done
}

#==========================================
# ESTRUCTURA DE CARPETAS
#==========================================
setup_directory_structure() {
    log_step "Creando estructura de carpetas"
    RUNDIR="$INITIAL_DIR/$WORKDIR/$(basename "${INPUT_PDB%.pdb}")_$(date +%Y%m%d_%H%M%S)"
    for d in 00_setup 01_minimization 02_equilibration 03_pulling 04_frames \
             05_umbrella_npt 06_umbrella_prod 07_analysis logs mdp_used; do
        create_dir "$RUNDIR/$d"
    done
    cp "$INITIAL_DIR/$BASE_MDP"/*.mdp "$RUNDIR/mdp_used/"

    # Guardar configuración
    cat > "$RUNDIR/config.txt" <<CONF
# Umbrella Sampling v2 — $(date)
SYSTEM_MODE=$SYSTEM_MODE
INPUT_PDB=$(basename "$INPUT_PDB")
CHAIN_A_RESIDUES=$CHAIN_A_RESIDUES
CHAIN_B_RESIDUES=$CHAIN_B_RESIDUES
REF_CHAIN_ITP=${REF_CHAIN_ITP:-}
PULL_GROUP1_NAME=${PULL_GROUP1_NAME:-Chain_A}
PULL_GROUP2_NAME=${PULL_GROUP2_NAME:-Chain_B}
BOX_CENTER=$BOX_CENTER
BOX_DIMS=$BOX_DIMS
PULL_DISTANCE=${PULL_DISTANCE:-5.0}
BOX_MARGIN=${BOX_MARGIN:-1.5}
PULL_DIM=$PULL_DIM
PULL_AXIS=${PULL_AXIS:-z}
PULL_RATE=$PULL_RATE
PULL_K=$PULL_K
PULL_TIME=$PULL_TIME
ION_CONC=$ION_CONC
WATER_MODEL=$WATER_MODEL
WINDOW_SPACING=$WINDOW_SPACING
FRAME_INCREMENT=$FRAME_INCREMENT
UMBRELLA_NS=$UMBRELLA_NS
PMF_UNIT=$PMF_UNIT
PARALLEL_WINDOWS=$PARALLEL_WINDOWS
CONF
    log_success "Directorio: $RUNDIR"
}

copy_common_files() {
    local t="$1"
    cp -f "$RUNDIR/00_setup/topol.top" "$t/"
    cp -f "$RUNDIR/00_setup/index.ndx" "$t/" 2>/dev/null || true
    for itp in "$RUNDIR/00_setup"/*.itp; do [ -f "$itp" ] && cp -f "$itp" "$t/"; done
    for ff in "$RUNDIR/00_setup"/*.ff; do [ -d "$ff" ] && cp -r "$ff" "$t/"; done
}

#==========================================
# ETAPA 1: TOPOLOGÍA
#==========================================
prepare_topology() {
    local STAGE="01_topology"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }

    log_step "Etapa 1: Preparando topología (pdb2gmx)"
    cd "$RUNDIR/00_setup" || exit 1
    cp "$INPUT_PDB" .
    local pdb_base; pdb_base=$(basename "$INPUT_PDB")

    if [ -n "$FF_NUMBER" ]; then
        echo "$FF_NUMBER" | run_gmx pdb2gmx -f "$pdb_base" -ignh -ter -o complex.gro \
            &> "$RUNDIR/logs/pdb2gmx.log" || {
            log_warning "Reintentando interactivamente..."
            run_gmx pdb2gmx -f "$pdb_base" -ignh -ter -o complex.gro 2>&1 | tee "$RUNDIR/logs/pdb2gmx.log"
        }
    else
        run_gmx pdb2gmx -f "$pdb_base" -ignh -ter -o complex.gro 2>&1 | tee "$RUNDIR/logs/pdb2gmx.log"
    fi
    [ ! -f complex.gro ] && { log_error "complex.gro no generado"; exit 1; }
    log_success "Topología generada"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 2: RESTRAINTS CADENA REFERENCIA
#==========================================
add_chain_restraints() {
    local STAGE="02_restraints"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }

    # Solo para modo multichain (proteína-ligando y membrana no requieren POSRES_B)
    if [ "$SYSTEM_MODE" != "multichain" ]; then
        log_info "Modo $SYSTEM_MODE: no requiere POSRES_B"
        mark_stage_done "$STAGE"; return
    fi

    log_step "Etapa 2: Position restraints cadena referencia"
    cd "$RUNDIR/00_setup" || exit 1

    local itp_file="topol_${REF_CHAIN_ITP}.itp"
    if [ ! -f "$itp_file" ]; then
        log_warning "$itp_file no encontrado, buscando..."
        ls -1 *.itp 2>/dev/null || true
        echo -e "  ${YELLOW}Nombre del .itp de la cadena de referencia:${NC}"
        read -r itp_file
        [ ! -f "$itp_file" ] && { log_error "$itp_file no encontrado"; exit 1; }
    fi

    local posre_file="posre_${REF_CHAIN_ITP}.itp"
    [ ! -f "$posre_file" ] && posre_file=$(ls posre_*.itp 2>/dev/null | grep -i "chain_b" | head -1 || true)
    [ -z "$posre_file" ] && posre_file=$(ls posre_*.itp 2>/dev/null | tail -1 || true)

    if [ -n "$posre_file" ] && ! grep -q "POSRES_B" "$itp_file"; then
        cat >> "$itp_file" <<EOF

; Position restraints for reference chain
#ifdef POSRES_B
#include "$posre_file"
#endif
EOF
        log_success "POSRES_B añadido a $itp_file"
    fi
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 3: DEFINIR CAJA (con auto-detección)
#==========================================
define_box() {
    local STAGE="03_box"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }

    log_step "Etapa 3: Definiendo caja de simulación"
    cd "$RUNDIR/00_setup" || exit 1

    if [ "$BOX_CENTER" = "auto" ]; then
        log_info "Auto-detectando dimensiones de caja..."

        # Obtener dimensiones del sistema
        local sys_info
        sys_info=$(run_gmx editconf -f complex.gro -o /dev/null -d 0 2>&1 || true)

        # Leer coordenadas extremas del GRO para calcular centro
        local xmin xmax ymin ymax zmin zmax
        xmin=$(awk 'NR>2 && NF>=6 {x=substr($0,21,8)+0; if(NR==3||x<min)min=x} END{print min}' complex.gro)
        xmax=$(awk 'NR>2 && NF>=6 {x=substr($0,21,8)+0; if(NR==3||x>max)max=x} END{print max}' complex.gro)
        ymin=$(awk 'NR>2 && NF>=6 {y=substr($0,29,8)+0; if(NR==3||y<min)min=y} END{print min}' complex.gro)
        ymax=$(awk 'NR>2 && NF>=6 {y=substr($0,29,8)+0; if(NR==3||y>max)max=y} END{print max}' complex.gro)
        zmin=$(awk 'NR>2 && NF>=6 {z=substr($0,37,8)+0; if(NR==3||z<min)min=z} END{print min}' complex.gro)
        zmax=$(awk 'NR>2 && NF>=6 {z=substr($0,37,8)+0; if(NR==3||z>max)max=z} END{print max}' complex.gro)

        # Centro del sistema
        local cx cy cz
        cx=$(echo "($xmin + $xmax) / 2" | bc -l)
        cy=$(echo "($ymin + $ymax) / 2" | bc -l)
        cz=$(echo "($zmin + $zmax) / 2" | bc -l)

        # Tamaño del sistema
        local sx sy sz
        sx=$(echo "$xmax - $xmin" | bc -l)
        sy=$(echo "$ymax - $ymin" | bc -l)
        sz=$(echo "$zmax - $zmin" | bc -l)

        # Dimensiones: sistema + 2*margen en ejes no-pulling
        # En eje de pulling: sistema + 2*(pull_distance + margen)
        local bx by bz m pd
        m=${BOX_MARGIN:-1.5}
        pd=${PULL_DISTANCE:-5.0}

        case $PULL_AXIS in
            x) bx=$(echo "$sx + 2*($pd + $m)" | bc -l)
               by=$(echo "$sy + 2*$m" | bc -l)
               bz=$(echo "$sz + 2*$m" | bc -l) ;;
            y) bx=$(echo "$sx + 2*$m" | bc -l)
               by=$(echo "$sy + 2*($pd + $m)" | bc -l)
               bz=$(echo "$sz + 2*$m" | bc -l) ;;
            *) bx=$(echo "$sx + 2*$m" | bc -l)
               by=$(echo "$sy + 2*$m" | bc -l)
               bz=$(echo "$sz + 2*($pd + $m)" | bc -l) ;;
        esac

        BOX_CENTER=$(printf "%.3f %.3f %.3f" "$cx" "$cy" "$cz")
        BOX_DIMS=$(printf "%.3f %.3f %.3f" "$bx" "$by" "$bz")
        log_success "Auto-caja: center=($BOX_CENTER) dims=($BOX_DIMS)"
    fi

    # shellcheck disable=SC2086
    run_gmx editconf -f complex.gro -o newbox.gro -center $BOX_CENTER -box $BOX_DIMS \
        &> "$RUNDIR/logs/editconf_box.log"
    [ ! -f newbox.gro ] && { log_error "newbox.gro no generado"; exit 1; }
    log_success "Caja definida: ($BOX_DIMS)"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPAS 4-6: SOLVATAR, MINIMIZAR, NPT
#==========================================
solvate_and_ionize() {
    local STAGE="04_solvate"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 4: Solvatación e iones"
    cd "$RUNDIR/00_setup" || exit 1
    run_gmx solvate -cp newbox.gro -cs "$WATER_FILE" -o solv.gro -p topol.top &> "$RUNDIR/logs/solvate.log"
    run_gmx grompp -f "$RUNDIR/mdp_used/ions.mdp" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 &> "$RUNDIR/logs/grompp_ions.log"
    echo "SOL" | run_gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
        -pname NA -nname CL -neutral -conc "$ION_CONC" &> "$RUNDIR/logs/genion.log"
    log_success "Iones añadidos ($ION_CONC M)"
    mark_stage_done "$STAGE"
}

run_minimization() {
    local STAGE="05_minimization"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 5: Minimización"
    cd "$RUNDIR/01_minimization" || exit 1; copy_common_files .
    cp "$RUNDIR/00_setup/solv_ions.gro" .
    run_gmx grompp -f "$RUNDIR/mdp_used/minim.mdp" -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2 &> "$RUNDIR/logs/grompp_em.log"
    $MDRUN -deffnm em -v &> "$RUNDIR/logs/mdrun_em.log"
    [ ! -f em.gro ] && { log_error "Minimización falló"; exit 1; }
    log_success "Minimización completada"
    mark_stage_done "$STAGE"
}

run_npt_equilibration() {
    local STAGE="06_npt"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 6: Equilibración NPT"
    cd "$RUNDIR/02_equilibration" || exit 1; copy_common_files .
    cp "$RUNDIR/01_minimization/em.gro" .
    run_gmx grompp -f "$RUNDIR/mdp_used/npt.mdp" -c em.gro -p topol.top -r em.gro -o npt.tpr -maxwarn 3 &> "$RUNDIR/logs/grompp_npt.log"
    $MDRUN -deffnm npt -v &> "$RUNDIR/logs/mdrun_npt.log"
    [ ! -f npt.gro ] && { log_error "NPT falló"; exit 1; }
    log_success "NPT completada"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 7: GRUPOS DE ÍNDICE (según modo)
#==========================================
create_pull_groups() {
    local STAGE="07_index"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 7: Grupos de índice ($SYSTEM_MODE)"
    cd "$RUNDIR/02_equilibration" || exit 1

    local last_idx
    last_idx=$(echo q | "$GMX" make_ndx -f npt.gro 2>&1 | grep -oP '^\s*\K\d+(?=\s)' | tail -1)
    local g1=$((last_idx + 1)); local g2=$((last_idx + 2))

    case $SYSTEM_MODE in
        multichain)
            run_gmx make_ndx -f npt.gro -o index.ndx &> "$RUNDIR/logs/make_ndx.log" <<EOF
r $CHAIN_A_RESIDUES
name $g1 Chain_A
r $CHAIN_B_RESIDUES
name $g2 Chain_B
q
EOF
            ;;
        protlig)
            run_gmx make_ndx -f npt.gro -o index.ndx &> "$RUNDIR/logs/make_ndx.log" <<EOF
r ${LIG_RESNAME:-LIG}
name $g1 LIG
1
name $g2 Protein_ref
q
EOF
            PULL_GROUP1_NAME="LIG"; PULL_GROUP2_NAME="Protein_ref"
            ;;
        membrane)
            run_gmx make_ndx -f npt.gro -o index.ndx &> "$RUNDIR/logs/make_ndx.log" <<EOF
r $CHAIN_A_RESIDUES
name $g1 Permeant
r ${MEMBRANE_LIPID:-POPC}
name $g2 Membrane
q
EOF
            PULL_GROUP1_NAME="Permeant"; PULL_GROUP2_NAME="Membrane"
            ;;
    esac

    cp index.ndx "$RUNDIR/00_setup/"
    log_success "Index creado: $PULL_GROUP1_NAME + $PULL_GROUP2_NAME"

    # Actualizar MDPs con los nombres de grupo correctos
    for mdp_file in "$RUNDIR/mdp_used"/md_pull.mdp "$RUNDIR/mdp_used"/npt_umbrella.mdp "$RUNDIR/mdp_used"/md_umbrella.mdp; do
        sed -i "s/^pull_group1_name.*/pull_group1_name        = $PULL_GROUP1_NAME/" "$mdp_file"
        sed -i "s/^pull_group2_name.*/pull_group2_name        = $PULL_GROUP2_NAME/" "$mdp_file"
    done
    log_success "MDPs actualizados con grupos: $PULL_GROUP1_NAME / $PULL_GROUP2_NAME"

    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 8: PULLING (SMD)
#==========================================
run_pulling() {
    local STAGE="08_pulling"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 8: Simulación de pulling (Steered MD)"
    cd "$RUNDIR/03_pulling" || exit 1; copy_common_files .
    cp "$RUNDIR/02_equilibration/npt.gro" .
    cp "$RUNDIR/02_equilibration/npt.cpt" .
    cp "$RUNDIR/02_equilibration/index.ndx" .

    local pull_mdp="$RUNDIR/mdp_used/md_pull.mdp"
    backup_mdp "$pull_mdp"
    sed -i "s/^nsteps.*/nsteps      = $PULL_NSTEPS    ; $PULL_TIME ps/" "$pull_mdp"
    sed -i "s/^pull_coord1_rate.*/pull_coord1_rate        = $PULL_RATE/" "$pull_mdp"
    sed -i "s/^pull_coord1_k.*/pull_coord1_k           = $PULL_K/" "$pull_mdp"
    sed -i "s/^pull_coord1_dim.*/pull_coord1_dim         = $PULL_DIM/" "$pull_mdp"

    run_gmx grompp -f "$pull_mdp" -c npt.gro -p topol.top -r npt.gro \
        -n index.ndx -t npt.cpt -o pull.tpr -maxwarn 3 &> "$RUNDIR/logs/grompp_pull.log"
    $MDRUN -deffnm pull -pf pullf.xvg -px pullx.xvg -v &> "$RUNDIR/logs/mdrun_pull.log"
    [ ! -f pull.xtc ] && { log_error "Pulling falló"; exit 1; }
    log_success "Pulling completado"
    mark_stage_done "$STAGE"
}

#==========================================
# CORRECCIÓN AUTOMÁTICA DE PBC
#==========================================
correct_pbc() {
    local STAGE="08b_pbc"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 8b: Corrección de PBC en trayectoria de pulling"
    cd "$RUNDIR/03_pulling" || exit 1

    # Centering: nojump + center + compact
    log_info "Aplicando corrección nojump..."
    echo "0" | run_gmx trjconv -s pull.tpr -f pull.xtc -o pull_nojump.xtc \
        -pbc nojump &> "$RUNDIR/logs/pbc_nojump.log" || {
        log_warning "PBC nojump falló, usando trayectoria original"
        mark_stage_done "$STAGE"; return
    }

    log_info "Centrando trayectoria..."
    echo -e "1\n0" | run_gmx trjconv -s pull.tpr -f pull_nojump.xtc -o pull_centered.xtc \
        -pbc mol -center &> "$RUNDIR/logs/pbc_center.log" || {
        log_warning "PBC center falló, usando nojump"
        cp pull_nojump.xtc pull_corrected.xtc
        mark_stage_done "$STAGE"; return
    }

    mv pull_centered.xtc pull_corrected.xtc
    rm -f pull_nojump.xtc

    log_success "PBC corregido → pull_corrected.xtc"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 9: EXTRACCIÓN DE FRAMES Y DISTANCIAS
#==========================================
extract_frames_and_distances() {
    local STAGE="09_frames"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 9: Extrayendo frames y calculando distancias COM"
    cd "$RUNDIR/04_frames" || exit 1
    cp "$RUNDIR/03_pulling/pull.tpr" .
    cp "$RUNDIR/02_equilibration/index.ndx" .

    # Usar trayectoria corregida si existe
    if [ -f "$RUNDIR/03_pulling/pull_corrected.xtc" ]; then
        cp "$RUNDIR/03_pulling/pull_corrected.xtc" pull.xtc
        log_info "Usando trayectoria con PBC corregido"
    else
        cp "$RUNDIR/03_pulling/pull.xtc" .
    fi

    echo 0 | run_gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep \
        &> "$RUNDIR/logs/trjconv_extract.log"
    local num_frames; num_frames=$(ls conf*.gro 2>/dev/null | wc -l)
    log_success "Frames extraídos: $num_frames"

    log_info "Calculando distancias COM..."
    : > summary_distances.dat
    for (( i=0; i<num_frames; i++ )); do
        if [ -f "conf${i}.gro" ]; then
            run_gmx distance -s pull.tpr -f "conf${i}.gro" -n index.ndx \
                -select "com of group \"$PULL_GROUP1_NAME\" plus com of group \"$PULL_GROUP2_NAME\"" \
                -oall "dist${i}.xvg" &> /dev/null || continue
            local d; d=$(tail -n 1 "dist${i}.xvg" 2>/dev/null | awk '{print $2}')
            [ -n "$d" ] && echo "${i} ${d}" >> summary_distances.dat
            rm -f "dist${i}.xvg"
        fi
        (( i % 50 == 0 )) && echo -ne "\r  Frame $i / $num_frames..."
    done
    echo ""
    [ ! -s summary_distances.dat ] && { log_error "No se calcularon distancias"; exit 1; }
    log_success "Distancias en summary_distances.dat"
    mark_stage_done "$STAGE"
}

#==========================================
# ESTIMACIÓN JARZYNSKI ΔG DESDE SMD
#==========================================
estimate_jarzynski() {
    local STAGE="09b_jarzynski"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 9b: Estimación Jarzynski ΔG desde SMD"

    local pullf="$RUNDIR/03_pulling/pullf.xvg"
    local pullx="$RUNDIR/03_pulling/pullx.xvg"

    if [ ! -f "$pullf" ] || [ ! -f "$pullx" ]; then
        log_warning "Archivos pullf/pullx no encontrados, saltando Jarzynski"
        mark_stage_done "$STAGE"; return
    fi

    python3 "$SCRIPT_DIR/jarzynski.py" "$pullf" "$pullx" "$RUNDIR/07_analysis" "$PMF_UNIT"

    log_success "Estimación Jarzynski guardada en 07_analysis/"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 10: SELECCIÓN ADAPTATIVA DE VENTANAS
#==========================================
select_windows() {
    local STAGE="10_windows"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 10: Selección adaptativa de ventanas"
    cd "$RUNDIR/04_frames" || exit 1

    log_info "Selección adaptativa (espaciado ≈ $WINDOW_SPACING nm)..."

    python3 "$SCRIPT_DIR/select_windows.py" "$WINDOW_SPACING" \
        summary_distances.dat "$RUNDIR/selected_windows.dat"

    # Leer resultados
    NUM_WINDOWS=$(grep -c '^[0-9]' "$RUNDIR/selected_windows.dat" || echo 0)
    local dmin dmax
    dmin=$(grep '^[0-9]' "$RUNDIR/selected_windows.dat" | head -1 | awk '{print $3}')
    dmax=$(grep '^[0-9]' "$RUNDIR/selected_windows.dat" | tail -1 | awk '{print $3}')

    log_success "Ventanas: $NUM_WINDOWS (rango: $dmin — $dmax nm)"

    # Mostrar tabla
    echo -e "\n${CYAN}  Ventana  Frame   Dist (nm)${NC}"
    echo    "  ──────── ─────── ─────────"
    while IFS= read -r line; do
        [[ "$line" =~ ^# ]] && continue
        printf "  %4s     %5s   %s\n" $(echo "$line" | awk '{print $1, $2, $3}')
    done < "$RUNDIR/selected_windows.dat"
    echo ""

    echo -e "  ${YELLOW}¿Continuar? (s/n) [s]:${NC}"
    read -r CONFIRM
    [[ "${CONFIRM,,}" == "n" ]] && { log_info "Modifica WINDOW_SPACING y reintenta"; exit 0; }

    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 11: NPT UMBRELLA POR VENTANA
#==========================================
run_umbrella_npt() {
    local STAGE="11_umbrella_npt"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 11: NPT equilibración por ventana"

    local npt_mdp="$RUNDIR/mdp_used/npt_umbrella.mdp"
    backup_mdp "$npt_mdp"
    sed -i "s/^pull_coord1_dim.*/pull_coord1_dim         = $PULL_DIM/" "$npt_mdp"

    local wins=()
    while IFS= read -r line; do
        [[ "$line" =~ ^# ]] && continue
        wins+=("$line")
    done < "$RUNDIR/selected_windows.dat"

    local total=${#wins[@]}; local done_count=0

    run_single_npt() {
        local w=$1 frame=$2
        local wdir="$RUNDIR/05_umbrella_npt/window_${w}"
        [ -f "$wdir/.done" ] && return 0
        mkdir -p "$wdir"; cd "$wdir"
        copy_common_files .
        cp "$RUNDIR/04_frames/conf${frame}.gro" conf.gro

        # Checkpoint granular: si existe .cpt de esta ventana, reanudar
        if [ -f "npt${w}.cpt" ] && [ -f "npt${w}.tpr" ]; then
            log_info "Reanudando NPT ventana $w desde checkpoint"
            $MDRUN -deffnm "npt${w}" -cpi "npt${w}.cpt" -v \
                &> "$RUNDIR/logs/mdrun_npt_w${w}_resume.log"
        else
            run_gmx grompp -f "$npt_mdp" -c conf.gro -p topol.top -r conf.gro \
                -n "$RUNDIR/00_setup/index.ndx" -o npt${w}.tpr -maxwarn 3 \
                &> "$RUNDIR/logs/grompp_npt_w${w}.log"
            $MDRUN -deffnm "npt${w}" -v &> "$RUNDIR/logs/mdrun_npt_w${w}.log"
        fi
        [ -f "npt${w}.gro" ] && touch .done
    }

    for entry in "${wins[@]}"; do
        local w frame dist
        read -r w frame dist <<< "$entry"

        if [ "$PARALLEL_WINDOWS" -le 1 ]; then
            echo -ne "\r  Ventana $((done_count+1))/$total (w=$w, frame=$frame)..."
            run_single_npt "$w" "$frame"
        else
            run_single_npt "$w" "$frame" &
            while [ "$(jobs -rp | wc -l)" -ge "$PARALLEL_WINDOWS" ]; do sleep 2; done
        fi
        done_count=$((done_count + 1))
    done
    wait
    echo ""
    log_success "NPT umbrella: $total ventanas"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 12: PRODUCCIÓN UMBRELLA
#==========================================
run_umbrella_production() {
    local STAGE="12_umbrella_prod"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 12: Producción umbrella ($UMBRELLA_NS ns/ventana)"

    local prod_mdp="$RUNDIR/mdp_used/md_umbrella.mdp"
    backup_mdp "$prod_mdp"
    sed -i "s/^nsteps.*/nsteps      = $UMBRELLA_NSTEPS    ; $UMBRELLA_NS ns/" "$prod_mdp"
    sed -i "s/^pull_coord1_dim.*/pull_coord1_dim         = $PULL_DIM/" "$prod_mdp"

    local wins=()
    while IFS= read -r line; do
        [[ "$line" =~ ^# ]] && continue
        wins+=("$line")
    done < "$RUNDIR/selected_windows.dat"

    local total=${#wins[@]}; local done_count=0

    run_single_prod() {
        local w=$1
        local nptdir="$RUNDIR/05_umbrella_npt/window_${w}"
        local pdir="$RUNDIR/06_umbrella_prod"
        [ -f "$pdir/umbrella${w}.gro" ] && return 0
        cd "$pdir"
        copy_common_files .

        # Checkpoint granular: si existe .cpt de esta ventana, reanudar
        if [ -f "umbrella${w}.cpt" ] && [ -f "umbrella${w}.tpr" ]; then
            log_info "Reanudando producción ventana $w desde checkpoint"
            $MDRUN -deffnm "umbrella${w}" -cpi "umbrella${w}.cpt" \
                -pf "umbrella${w}_pullf.xvg" -px "umbrella${w}_pullx.xvg" -v \
                -append \
                &> "$RUNDIR/logs/mdrun_umb_w${w}_resume.log"
        else
            run_gmx grompp -f "$prod_mdp" -c "$nptdir/npt${w}.gro" \
                -p topol.top -r "$nptdir/npt${w}.gro" -t "$nptdir/npt${w}.cpt" \
                -n "$RUNDIR/00_setup/index.ndx" -o "umbrella${w}.tpr" -maxwarn 3 \
                &> "$RUNDIR/logs/grompp_umb_w${w}.log"
            $MDRUN -deffnm "umbrella${w}" -pf "umbrella${w}_pullf.xvg" \
                -px "umbrella${w}_pullx.xvg" -v \
                &> "$RUNDIR/logs/mdrun_umb_w${w}.log"
        fi
    }

    for entry in "${wins[@]}"; do
        local w frame dist
        read -r w frame dist <<< "$entry"
        if [ "$PARALLEL_WINDOWS" -le 1 ]; then
            echo -ne "\r  Ventana $((done_count+1))/$total (w=$w)..."
            run_single_prod "$w"
        else
            run_single_prod "$w" &
            while [ "$(jobs -rp | wc -l)" -ge "$PARALLEL_WINDOWS" ]; do sleep 2; done
        fi
        done_count=$((done_count + 1))
    done
    wait
    echo ""
    log_success "Producción umbrella: $total ventanas"
    mark_stage_done "$STAGE"
}

#==========================================
# ETAPA 13: WHAM + ANÁLISIS
#==========================================
run_wham_analysis() {
    local STAGE="13_wham"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 13: Análisis WHAM → PMF"
    cd "$RUNDIR/07_analysis" || exit 1

    ls "$RUNDIR/06_umbrella_prod"/umbrella*.tpr 2>/dev/null | sort -V > tpr-files.dat
    ls "$RUNDIR/06_umbrella_prod"/umbrella*_pullf.xvg 2>/dev/null | sort -V > pullf-files.dat

    local n_tpr; n_tpr=$(wc -l < tpr-files.dat)
    local n_pullf; n_pullf=$(wc -l < pullf-files.dat)

    [ "$n_tpr" -eq 0 ] || [ "$n_pullf" -eq 0 ] && { log_error "No hay archivos para WHAM"; exit 1; }

    run_gmx wham -it tpr-files.dat -if pullf-files.dat \
        -o pmf.xvg -hist histogram.xvg -unit "$PMF_UNIT" \
        -nBootstrap 200 -bs-method traj \
        -bsres bootstrap_results.xvg -bsprof bootstrap_profiles.xvg \
        &> "$RUNDIR/logs/wham.log"

    [ ! -f pmf.xvg ] && { log_error "WHAM falló"; exit 1; }

    cp "$RUNDIR/03_pulling/pullf.xvg" pull_force.xvg 2>/dev/null || true
    cp "$RUNDIR/04_frames/summary_distances.dat" . 2>/dev/null || true
    cp "$RUNDIR/selected_windows.dat" . 2>/dev/null || true

    log_success "WHAM completado: pmf.xvg, histogram.xvg"

    # Extraer ΔG
    local min_e max_e delta_g
    min_e=$(grep -v '^[@#]' pmf.xvg | awk '{print $2}' | sort -n | head -1)
    max_e=$(grep -v '^[@#]' pmf.xvg | awk '{print $2}' | sort -n | tail -1)
    delta_g=$(echo "$max_e - $min_e" | bc -l 2>/dev/null || echo "N/A")
    echo -e "\n  ${MAGENTA}══════════════════════════════${NC}"
    echo -e "  ${MAGENTA}  ΔG ≈ ${delta_g} ${PMF_UNIT}${NC}"
    echo -e "  ${MAGENTA}══════════════════════════════${NC}\n"

    mark_stage_done "$STAGE"
}

#==========================================
# ANÁLISIS DE CONVERGENCIA POR VENTANA
#==========================================
analyze_convergence() {
    local STAGE="14_convergence"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 14: Análisis de convergencia (block averaging)"
    cd "$RUNDIR/07_analysis" || exit 1

    python3 "$SCRIPT_DIR/analyze_convergence.py" \
        "$RUNDIR/06_umbrella_prod" "$RUNDIR/07_analysis" "$PMF_UNIT"

    log_success "Reporte en: convergence_report.dat"
    mark_stage_done "$STAGE"
}

#==========================================
# DETECCIÓN DE GAPS + AUTO-FILL
#==========================================
detect_histogram_gaps() {
    local STAGE="15_gaps"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 15: Detección de gaps en histogramas"
    cd "$RUNDIR/07_analysis" || exit 1

    [ ! -f histogram.xvg ] && { log_warning "histogram.xvg no encontrado"; mark_stage_done "$STAGE"; return; }

    python3 "$SCRIPT_DIR/detect_gaps.py" histogram.xvg "$WINDOW_SPACING" \
        "$RUNDIR/04_frames/summary_distances.dat" "$RUNDIR/gap_windows.dat"

    # Si hay gaps, preguntar si agregar ventanas
    if [ -f "$RUNDIR/gap_windows.dat" ] && grep -q '^[0-9]' "$RUNDIR/gap_windows.dat"; then
        local n_gaps; n_gaps=$(grep -c '^[0-9]' "$RUNDIR/gap_windows.dat")
        echo -e "\n  ${YELLOW}¿Ejecutar $n_gaps ventanas adicionales para llenar gaps? (s/n) [s]:${NC}"
        read -r FILL_GAPS
        if [[ "${FILL_GAPS,,}" != "n" ]]; then
            log_info "Ejecutando ventanas gap-fill..."
            run_gap_fill_windows
        fi
    fi

    mark_stage_done "$STAGE"
}

run_gap_fill_windows() {
    local npt_mdp="$RUNDIR/mdp_used/npt_umbrella.mdp"
    local prod_mdp="$RUNDIR/mdp_used/md_umbrella.mdp"
    local gap_idx=0

    while IFS= read -r line; do
        [[ "$line" =~ ^# ]] && continue
        local frame dist
        read -r frame dist _ <<< "$line"

        local gw="gap${gap_idx}"
        log_info "Gap-fill $gw: frame=$frame, dist=$dist"

        # NPT
        local wdir="$RUNDIR/05_umbrella_npt/window_${gw}"
        mkdir -p "$wdir"; cd "$wdir"
        copy_common_files .
        cp "$RUNDIR/04_frames/conf${frame}.gro" conf.gro
        run_gmx grompp -f "$npt_mdp" -c conf.gro -p topol.top -r conf.gro \
            -n "$RUNDIR/00_setup/index.ndx" -o "npt${gw}.tpr" -maxwarn 3 \
            &> "$RUNDIR/logs/grompp_npt_${gw}.log"
        $MDRUN -deffnm "npt${gw}" -v &> "$RUNDIR/logs/mdrun_npt_${gw}.log"

        # Producción
        cd "$RUNDIR/06_umbrella_prod"
        copy_common_files .
        run_gmx grompp -f "$prod_mdp" -c "$wdir/npt${gw}.gro" \
            -p topol.top -r "$wdir/npt${gw}.gro" -t "$wdir/npt${gw}.cpt" \
            -n "$RUNDIR/00_setup/index.ndx" -o "umbrella${gw}.tpr" -maxwarn 3 \
            &> "$RUNDIR/logs/grompp_umb_${gw}.log"
        $MDRUN -deffnm "umbrella${gw}" -pf "umbrella${gw}_pullf.xvg" \
            -px "umbrella${gw}_pullx.xvg" -v \
            &> "$RUNDIR/logs/mdrun_umb_${gw}.log"

        gap_idx=$((gap_idx + 1))
    done < "$RUNDIR/gap_windows.dat"

    # Re-ejecutar WHAM con ventanas adicionales
    log_info "Re-ejecutando WHAM con ventanas gap-fill..."
    cd "$RUNDIR/07_analysis"
    ls "$RUNDIR/06_umbrella_prod"/umbrella*.tpr 2>/dev/null | sort -V > tpr-files.dat
    ls "$RUNDIR/06_umbrella_prod"/umbrella*_pullf.xvg 2>/dev/null | sort -V > pullf-files.dat

    run_gmx wham -it tpr-files.dat -if pullf-files.dat \
        -o pmf.xvg -hist histogram.xvg -unit "$PMF_UNIT" \
        -nBootstrap 200 -bs-method traj \
        -bsres bootstrap_results.xvg -bsprof bootstrap_profiles.xvg \
        &> "$RUNDIR/logs/wham_gapfill.log"

    log_success "WHAM re-ejecutado con ventanas gap-fill"
}

#==========================================
# VALIDACIÓN DE OVERLAP
#==========================================
validate_overlap() {
    local STAGE="15b_overlap"
    is_stage_done "$STAGE" && { log_info "Saltando: $STAGE"; return; }
    log_step "Etapa 15b: Validación cuantitativa de overlap"
    cd "$RUNDIR/07_analysis" || exit 1

    [ ! -f histogram.xvg ] && { log_warning "histogram.xvg no encontrado"; mark_stage_done "$STAGE"; return; }

    python3 "$SCRIPT_DIR/validate_overlap.py" histogram.xvg "$RUNDIR/07_analysis"

    log_success "Reporte en: overlap_report.dat"
    mark_stage_done "$STAGE"
}

#==========================================
# LIMPIEZA DE TEMPORALES
#==========================================
cleanup_intermediates() {
    local target_dir="${1:-$RUNDIR}"
    log_step "Limpieza de archivos temporales"

    if [ ! -d "$target_dir" ]; then
        log_error "Directorio no encontrado: $target_dir"; return 1
    fi

    # Calcular espacio ocupado antes
    local size_before
    size_before=$(du -sh "$target_dir" 2>/dev/null | cut -f1)
    log_info "Tamaño actual: $size_before"

    # Contar archivos temporales
    local n_conf n_dist n_nojump n_backup
    n_conf=$(find "$target_dir/04_frames" -name 'conf*.gro' 2>/dev/null | wc -l)
    n_dist=$(find "$target_dir/04_frames" -name 'dist*.xvg' 2>/dev/null | wc -l)
    n_nojump=$(find "$target_dir/03_pulling" -name 'pull_nojump.xtc' 2>/dev/null | wc -l)
    n_backup=$(find "$target_dir" -name '*.mdp.original' 2>/dev/null | wc -l)

    echo -e "\n  ${CYAN}Archivos temporales encontrados:${NC}"
    echo "    conf*.gro (frames extraídos):  $n_conf archivos"
    echo "    dist*.xvg (distancias temp):   $n_dist archivos"
    echo "    pull_nojump.xtc (PBC inter):   $n_nojump archivos"
    echo "    *.mdp.original (backups MDP):  $n_backup archivos"
    echo ""

    local total=$((n_conf + n_dist + n_nojump))
    if [ "$total" -eq 0 ]; then
        log_info "No hay archivos temporales para limpiar"
        return 0
    fi

    echo -e "  ${CYAN}Opciones de limpieza:${NC}"
    echo "    1) Comprimir frames (conf*.gro → frames.tar.gz)"  
    echo "    2) Eliminar frames y temporales (conservar summary_distances.dat)"
    echo "    3) Eliminar todo lo temporal + comprimir logs"
    echo "    0) Cancelar"
    echo -e "  ${YELLOW}Seleccione [0-3] (default: 1):${NC}"
    read -r CLEAN_MODE
    CLEAN_MODE=${CLEAN_MODE:-1}

    case $CLEAN_MODE in
        1)
            log_info "Comprimiendo frames..."
            if [ "$n_conf" -gt 0 ]; then
                tar -czf "$target_dir/04_frames/frames_backup.tar.gz" \
                    -C "$target_dir/04_frames" $(cd "$target_dir/04_frames" && ls conf*.gro 2>/dev/null) \
                    2>/dev/null && {
                    rm -f "$target_dir/04_frames"/conf*.gro
                    log_success "Frames comprimidos → frames_backup.tar.gz"
                }
            fi
            rm -f "$target_dir/04_frames"/dist*.xvg 2>/dev/null
            rm -f "$target_dir/03_pulling/pull_nojump.xtc" 2>/dev/null
            ;;
        2)
            log_info "Eliminando archivos temporales..."
            rm -f "$target_dir/04_frames"/conf*.gro 2>/dev/null
            rm -f "$target_dir/04_frames"/dist*.xvg 2>/dev/null
            rm -f "$target_dir/03_pulling/pull_nojump.xtc" 2>/dev/null
            log_success "Temporales eliminados"
            ;;
        3)
            log_info "Limpieza profunda..."
            rm -f "$target_dir/04_frames"/conf*.gro 2>/dev/null
            rm -f "$target_dir/04_frames"/dist*.xvg 2>/dev/null
            rm -f "$target_dir/03_pulling/pull_nojump.xtc" 2>/dev/null
            # Comprimir logs
            if [ -d "$target_dir/logs" ]; then
                tar -czf "$target_dir/logs_backup.tar.gz" \
                    -C "$target_dir" logs/ 2>/dev/null && {
                    rm -rf "$target_dir/logs"
                    log_success "Logs comprimidos → logs_backup.tar.gz"
                }
            fi
            # Comprimir trayectorias de umbrella (guardar solo .tpr y pullf/pullx)
            for xtc in "$target_dir/06_umbrella_prod"/umbrella*.xtc; do
                [ -f "$xtc" ] && rm -f "$xtc"
            done
            log_success "Limpieza profunda completada"
            ;;
        0) log_info "Limpieza cancelada"; return 0 ;;
        *) log_warning "Opción inválida"; return 0 ;;
    esac

    local size_after
    size_after=$(du -sh "$target_dir" 2>/dev/null | cut -f1)
    echo -e "\n  ${GREEN}Espacio: $size_before → $size_after${NC}"
}

#==========================================
# RESUMEN FINAL
#==========================================
generate_summary() {
    log_step "Generando resumen"
    cat > "$RUNDIR/RESUMEN.txt" <<SUMMARY
═══════════════════════════════════════
  UMBRELLA SAMPLING v3 — RESUMEN
  $(date)
═══════════════════════════════════════

Modo:           $SYSTEM_MODE
Estructura:     $(basename "$INPUT_PDB")
Grupo móvil:    $CHAIN_A_RESIDUES ($PULL_GROUP1_NAME)
Grupo ref:      $CHAIN_B_RESIDUES ($PULL_GROUP2_NAME)

Caja: center=($BOX_CENTER) dims=($BOX_DIMS)
Agua: $WATER_MODEL | Iones: $ION_CONC M

Pulling: rate=$PULL_RATE nm/ps, k=$PULL_K kJ/mol/nm², t=$PULL_TIME ps
Dirección: $PULL_DIM

Umbrella:
  Ventanas: $(grep -c '^[0-9]' "$RUNDIR/selected_windows.dat" 2>/dev/null || echo "N/A")
  Espaciado: $WINDOW_SPACING nm (adaptativo)
  Producción: $UMBRELLA_NS ns/ventana
  PMF: $PMF_UNIT

Resultados: $RUNDIR/07_analysis/
SUMMARY
    log_success "Resumen en RESUMEN.txt"
}

#==========================================
# MAIN
#==========================================
main() {
    [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ] && show_help

    # Verificar que scripts/ existe
    if [ ! -d "$SCRIPT_DIR" ]; then
        log_error "Carpeta scripts/ no encontrada en: $SCRIPT_DIR"
        log_info "Asegúrate de que la carpeta scripts/ está junto a umbrella_pipeline.sh"
        exit 1
    fi

    # Modo --cleanup
    if [ "${1:-}" = "--cleanup" ]; then
        [ -z "${2:-}" ] || [ ! -d "${2:-}" ] && { log_error "Directorio inválido"; exit 1; }
        cleanup_intermediates "$(realpath "$2")"
        trap - EXIT
        exit 0
    fi

    if [ "${1:-}" = "--resume" ]; then
        [ -z "${2:-}" ] || [ ! -d "${2:-}" ] && { log_error "Directorio inválido"; exit 1; }
        RUNDIR="$(realpath "$2")"
        log_info "Reanudando desde: $RUNDIR"
        if [ -f "$RUNDIR/config.txt" ]; then
            # shellcheck disable=SC1090
            source <(grep -v '^#' "$RUNDIR/config.txt" | grep '=')
            INPUT_PDB="$INITIAL_DIR/$INPUT_PDB"
            PULL_NSTEPS=$(echo "$PULL_TIME / 0.002" | bc)
            UMBRELLA_NSTEPS=$(echo "$UMBRELLA_NS * 500000" | bc)
            # Defaults para variables v3
            PULL_GROUP1_NAME=${PULL_GROUP1_NAME:-Chain_A}
            PULL_GROUP2_NAME=${PULL_GROUP2_NAME:-Chain_B}
            SYSTEM_MODE=${SYSTEM_MODE:-multichain}
            PULL_AXIS=${PULL_AXIS:-z}
            log_success "Configuración cargada"
        else
            log_error "config.txt no encontrado"; exit 1
        fi
    else
        get_user_input
        validate_mdp_files
        setup_directory_structure
    fi

    # Pipeline principal
    prepare_topology
    add_chain_restraints
    define_box
    solvate_and_ionize
    run_minimization
    run_npt_equilibration
    create_pull_groups
    run_pulling
    correct_pbc                    # Corrección PBC
    extract_frames_and_distances
    estimate_jarzynski             # Jarzynski ΔG
    select_windows                 # Selección adaptativa
    run_umbrella_npt
    run_umbrella_production
    run_wham_analysis
    analyze_convergence            # Convergencia
    validate_overlap               # NUEVO: validación de overlap
    detect_histogram_gaps          # Detección + auto-fill
    generate_summary

    echo ""
    echo -e "${GREEN}═══════════════════════════════════════════${NC}"
    echo -e "${GREEN}  ✓ UMBRELLA SAMPLING v3 COMPLETADO${NC}"
    echo -e "${GREEN}═══════════════════════════════════════════${NC}"
    echo ""
    echo -e "  PMF:          ${CYAN}$RUNDIR/07_analysis/pmf.xvg${NC}"
    echo -e "  Histogramas:  ${CYAN}$RUNDIR/07_analysis/histogram.xvg${NC}"
    echo -e "  Overlap:      ${CYAN}$RUNDIR/07_analysis/overlap_report.dat${NC}"
    echo -e "  Convergencia: ${CYAN}$RUNDIR/07_analysis/convergence_report.dat${NC}"
    echo -e "  Jarzynski:    ${CYAN}$RUNDIR/07_analysis/jarzynski.dat${NC}"
    echo ""
    echo -e "  Graficar: ${YELLOW}python3 plot_umbrella.py $RUNDIR/07_analysis/${NC}"
    echo -e "  Limpiar:  ${YELLOW}$0 --cleanup $RUNDIR${NC}"
    echo ""

    # Ofrecer limpieza
    echo -e "  ${YELLOW}¿Limpiar archivos temporales ahora? (s/n) [n]:${NC}"
    read -r DO_CLEANUP
    if [[ "${DO_CLEANUP,,}" == "s" ]]; then
        cleanup_intermediates "$RUNDIR"
    fi

    trap - EXIT
}

main "$@"
