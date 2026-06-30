#!/usr/bin/env bash
set -euo pipefail

# Script para preparar el directorio de muestras y análisis

DATA_DIR="$HOME/lrs_benchmarking/samples_try"
ANALYSIS_DIR="$DATA_DIR/alignment_analysis"

# Crear estructura de directorios de análisis
mkdir -p "$ANALYSIS_DIR"/{scripts,tables,figures,logs}

# Verificar si DATA_DIR existe
if [[ -d "$DATA_DIR" ]]; then
  echo "Data directory found: $DATA_DIR"
else
  echo "Data directory NOT found: $DATA_DIR" >&2
  exit 1
fi

echo "Analysis directory: $ANALYSIS_DIR"
