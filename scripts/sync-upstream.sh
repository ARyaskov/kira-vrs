#!/usr/bin/env bash
# Re-vendor the GA4GH VRS specification test material at an exact revision.
#
# Usage: scripts/sync-upstream.sh [<git ref>]      (default: the tag recorded below)
#
# Copies into crates/kira-vrs-validation/upstream/:
#   validation/{models,functions}.json   (regenerated from the YAML sources)
#   examples/*.json + test_definitions.yaml
#   schema/vrs/*  schema/gkm-core/*      (generated JSON Schemas)
#   LICENSE-APACHE-2.0, REVISION.md
#
# After running it, update the revision constants in crates/kira-vrs/src/spec.rs and
# crates/kira-vrs/Cargo.toml ([package.metadata.vrs]) and review the diff — model, digest and
# normalization changes upstream require code changes here.
set -euo pipefail

REF="${1:-2.1.0}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$ROOT/crates/kira-vrs-validation/upstream"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

git clone --quiet https://github.com/ga4gh/vrs.git "$WORK/vrs"
git -C "$WORK/vrs" checkout --quiet "$REF"
git -C "$WORK/vrs" submodule update --init --quiet
VRS_REV="$(git -C "$WORK/vrs" rev-parse HEAD)"
VRS_DATE="$(git -C "$WORK/vrs" log -1 --format=%cs HEAD)"
CORE_REV="$(git -C "$WORK/vrs/submodules/gkm-core" rev-parse HEAD)"
CORE_VER="$(grep -oE 'gkm-core/[0-9][^/]*' "$WORK/vrs/submodules/gkm-core/schema/gkm-core/gkm-core-source.yaml" | head -1 | cut -d/ -f2)"

rm -rf "$DEST"
mkdir -p "$DEST/validation" "$DEST/examples" "$DEST/schema/vrs" "$DEST/schema/gkm-core"

python3 - "$WORK/vrs" "$DEST" <<'EOF'
import json, sys, pathlib
import yaml  # pip install pyyaml
src, dest = map(pathlib.Path, sys.argv[1:3])
for name in ("models", "functions"):
    data = yaml.safe_load((src / "validation" / f"{name}.yaml").read_text())
    (dest / "validation" / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n")
EOF

cp "$WORK"/vrs/examples/json/*.json "$DEST/examples/"
cp "$WORK/vrs/tests/test_definitions.yaml" "$DEST/examples/"
cp "$WORK"/vrs/schema/vrs/json/* "$DEST/schema/vrs/"
cp "$WORK"/vrs/submodules/gkm-core/schema/gkm-core/json/* "$DEST/schema/gkm-core/"
cp "$WORK/vrs/LICENSE" "$DEST/LICENSE-APACHE-2.0"

cat > "$DEST/REVISION.md" <<EOF
# Vendored upstream revision

| Repository | Ref | Commit | Date |
|---|---|---|---|
| https://github.com/ga4gh/vrs | \`$REF\` | \`$VRS_REV\` | $VRS_DATE |
| https://github.com/ga4gh/gkm-core (submodule) | \`$CORE_VER\` | \`$CORE_REV\` | — |

Files under this directory are copied verbatim from those revisions (the validation JSON is
regenerated from the YAML sources) and are licensed under the Apache License 2.0
(see \`LICENSE-APACHE-2.0\`). Regenerate with \`scripts/sync-upstream.sh $REF\`.
EOF

echo "Vendored VRS $REF ($VRS_REV) with gkm-core $CORE_VER ($CORE_REV) into $DEST"
