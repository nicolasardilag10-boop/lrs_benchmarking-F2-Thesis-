#!/usr/bin/env bash

set -euo pipefail

VACMAP_REPO="https://github.com/micahvista/VACmap.git"
VACMAP_COMMIT="4a804758e3d54becdd81ef02c09ef8326f7e5830"

echo "=== VACMAP POST-DEPLOY ==="
echo "Conda prefix : $CONDA_PREFIX"
echo "Repository   : $VACMAP_REPO"
echo "Commit       : $VACMAP_COMMIT"

TMPDIR_VACMAP="$(mktemp -d)"

cleanup() {
    rm -rf "$TMPDIR_VACMAP"
}

trap cleanup EXIT

echo
echo "=== CLONING VACMAP ==="

git clone --quiet \
    "$VACMAP_REPO" \
    "$TMPDIR_VACMAP/VACmap"

git -C "$TMPDIR_VACMAP/VACmap" \
    checkout --quiet "$VACMAP_COMMIT"

ACTUAL_COMMIT="$(
    git -C "$TMPDIR_VACMAP/VACmap" rev-parse HEAD
)"

echo "Checked-out commit: $ACTUAL_COMMIT"

if [[ "$ACTUAL_COMMIT" != "$VACMAP_COMMIT" ]]; then
    echo "ERROR: VACmap commit mismatch."
    exit 1
fi

echo
echo "=== INSTALLING VACMAP ==="

"$CONDA_PREFIX/bin/python" \
    -m pip install \
    --no-build-isolation \
    --no-deps \
    "$TMPDIR_VACMAP/VACmap"

VACMAP_BIN="$CONDA_PREFIX/bin/vacmap"

if [[ ! -f "$VACMAP_BIN" ]]; then
    echo "ERROR: VACmap launcher was not installed."
    exit 1
fi

echo
echo "=== FIXING VACMAP LAUNCHER ==="

FIRST_LINE="$(head -n 1 "$VACMAP_BIN")"

if [[ "$FIRST_LINE" != '#!'* ]]; then
    TEMP_LAUNCHER="${VACMAP_BIN}.fixed"

    {
        printf '%s\n' '#!/usr/bin/env python3'
        cat "$VACMAP_BIN"
    } > "$TEMP_LAUNCHER"

    mv "$TEMP_LAUNCHER" "$VACMAP_BIN"
    chmod +x "$VACMAP_BIN"
fi

echo
echo "=== VACMAP VALIDATION ==="

"$VACMAP_BIN" --help >/dev/null

echo "VACMAP POST-DEPLOY: PASS"
