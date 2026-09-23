# Workflow namespace

This directory records the target organization for future workflow consolidation. Active workflows have not been moved here when doing so would break constitutional naming rules or relative paths.

Current entry points:

- aligners: root-level `ont.read_mapping.*.smk` and `pb.read_mapping.*.smk`;
- assemblers: [`../assemblers/`](../assemblers/README.md);
- variant callers: root-level caller `.smk` files and [`../SV aligners call/`](<../SV aligners call/README.md>).

Treat `workflows/` as a namespace and roadmap, not as the only executable source tree.
