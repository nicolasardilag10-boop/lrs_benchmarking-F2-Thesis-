#!/usr/bin/env bash

# Compatibility shim.
# Canonical shortcuts now live at:
#   ~/lrs_benchmarking/docs/lrs_shortcuts.sh

_LRS_COMPAT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${_LRS_COMPAT_DIR}/../docs/lrs_shortcuts.sh"
unset _LRS_COMPAT_DIR
