#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import re
import subprocess
from pathlib import Path, PurePosixPath

START = "<!-- AUTO_REPOSITORY_TREE_START -->"
END = "<!-- AUTO_REPOSITORY_TREE_END -->"

MODULE_READMES = {
    "alignment_analysis/README.md": "alignment_analysis",
    "assemblers/README.md": "assemblers",
    "assemblers/whole_genome_asm/README.md": "assemblers/whole_genome_asm",
    "assemblers/whole_genome_asm/assessment/README.md": "assemblers/whole_genome_asm/assessment",
    "assembly_analysis/README.md": "assembly_analysis",
    "variant_calling_analysis/README.md": "variant_calling_analysis",
}

ROOT_AREAS = [
    ("Alignment benchmark", "alignment_analysis", "metrics, canonical tables, and final figures"),
    ("Assembly workflows", "assemblers", "whole-genome and chromosome-21 assembly workflows"),
    ("Assembly analysis", "assembly_analysis", "assembly metric aggregation and figures"),
    ("Variant analysis", "variant_calling_analysis", "variant-calling analysis layer"),
    ("SV workspace", "SV aligners call", "SV and aligner validation workflows"),
    ("Final report", "final_report_files", "LaTeX report and integrated aligner workflow"),
    ("Documentation", "docs", "maps, methods, inventories, and troubleshooting"),
    ("Project figures", "figures", "cross-project and variant figures"),
    ("Results", "results", "organized links to generated outputs"),
    ("Tests", "tests", "fixtures, smoke tests, and validation"),
    ("Archive", "archive", "historical code, snapshots, and retired material"),
]

# Direct files in a top-level directory are listed only when there are few of them.
# Large result-heavy folders stay compact in the root README.
ROOT_DIRECT_FILE_LIMIT = 8

# These module directories start collapsed because they tend to be large.
COLLAPSED_BY_DEFAULT = {
    "figures",
    "tables",
    "results",
    "archive",
    "assessment",
    "scripts",
    "containers",
    "config",
    "docs",
}

SECTION_RE = re.compile(
    r"(?ms)^## Repository (?:structure|explorer)\n\n"
    r"<!-- AUTO_REPOSITORY_TREE_START -->.*?"
    r"^<!-- AUTO_REPOSITORY_TREE_END -->\n*"
)

def repo_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            text=True,
        ).strip()
    )

def repository_files(root: Path) -> list[str]:
    """Return visible version-control candidates that still exist on disk."""
    out = subprocess.check_output(
        [
            "git",
            "-C",
            str(root),
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
        ],
        text=True,
    )
    excluded_parts = {".git", ".snakemake", "__pycache__", "node_modules"}
    excluded_roots = {"overleaf-toolkit"}
    files = []

    for raw_path in out.splitlines():
        path = raw_path.strip()
        if not path:
            continue
        parts = PurePosixPath(path).parts
        if parts[0] in excluded_roots or excluded_parts.intersection(parts):
            continue
        if not (root / path).is_file():
            continue
        files.append(path)

    return sorted(set(files))

def esc(value: str) -> str:
    return html.escape(value, quote=False)

def md_file_link(target: str, label: str | None = None) -> str:
    label = label or PurePosixPath(target).name
    href = f"<{target}>" if " " in target else target
    return f"[`{esc(label)}`]({href})"

def immediate_children(files: list[str], prefix: str = ""):
    prefix = prefix.strip("/")
    needle = prefix + "/" if prefix else ""
    dirs = set()
    direct_files = []

    for path in files:
        if prefix and not path.startswith(needle):
            continue

        rel = path[len(needle):] if prefix else path
        if not rel:
            continue

        parts = rel.split("/")
        if len(parts) == 1:
            direct_files.append(parts[0])
        else:
            dirs.add(parts[0])

    return sorted(dirs), sorted(direct_files)

def count_under(files: list[str], prefix: str) -> int:
    needle = prefix.rstrip("/") + "/"
    return sum(1 for path in files if path.startswith(needle))

def render_root(files: list[str]) -> str:
    lines = [
        "## Repository explorer",
        "",
        START,
        "Generated from version-control candidates. This view highlights "
        "scientific entry points instead of listing every result file. See "
        "[`docs/REPOSITORY_TREE.md`](docs/REPOSITORY_TREE.md) for the expanded map.",
        "",
        "| Area | Location | Purpose |",
        "|---|---|---|",
    ]

    for label, dirname, purpose in ROOT_AREAS:
        readme = f"{dirname}/README.md"
        lines.append(
            f"| {esc(label)} | {md_file_link(readme, dirname + '/')} | {esc(purpose)} |"
        )

    root_files = {path for path in files if "/" not in path}
    workflow_groups = [
        ("Alignment", sorted(path for path in root_files if ".read_mapping." in path)),
        (
            "Small variants",
            sorted(path for path in root_files if ".snv_" in path or path == "run_happy.smk"),
        ),
        (
            "Structural variants",
            sorted(
                path
                for path in root_files
                if path
                in {
                    "cuteSV.hg38.smk",
                    "pbsv.hg38.smk",
                    "sawfish.hg38.smk",
                    "sniffles2.hg38.smk",
                    "truvari_anno.smk",
                }
            ),
        ),
    ]
    lines += ["", "<details>", "<summary><b>Constitution-required root workflows</b></summary>", ""]
    for label, paths in workflow_groups:
        if paths:
            links = ", ".join(md_file_link(path) for path in paths)
            lines.append(f"- **{label}:** {links}")
    lines += ["", "</details>", ""]

    lines.append(END)
    return "\n".join(lines)

def build_tree(files: list[str], module_root: str):
    prefix = module_root.rstrip("/") + "/"
    tree = {"dirs": {}, "files": []}

    for path in files:
        if not path.startswith(prefix):
            continue

        rel = path[len(prefix):]
        if not rel:
            continue

        node = tree
        parts = rel.split("/")

        for part in parts[:-1]:
            node = node["dirs"].setdefault(
                part,
                {"dirs": {}, "files": []},
            )

        node["files"].append(parts[-1])

    return tree

def recursive_count(node) -> int:
    return len(node["files"]) + sum(
        recursive_count(child)
        for child in node["dirs"].values()
    )

def render_node(node, rel_dir: str = "") -> list[str]:
    lines = []

    for filename in sorted(node["files"]):
        target = f"{rel_dir}/{filename}" if rel_dir else filename
        lines.append(f"- {md_file_link(target)}")

    for dirname in sorted(node["dirs"]):
        child = node["dirs"][dirname]
        child_rel = f"{rel_dir}/{dirname}" if rel_dir else dirname
        n = recursive_count(child)
        plural = "file" if n == 1 else "files"
        open_attr = "" if dirname in COLLAPSED_BY_DEFAULT else " open"

        lines += [
            "",
            f"<details{open_attr}>",
            f"<summary><b>{esc(dirname)}/</b> — {n} {plural}</summary>",
            "",
        ]

        lines.extend(render_node(child, child_rel) or ["_No tracked files._"])
        lines += ["", "</details>"]

    return lines

def render_module(files: list[str], module_root: str) -> str:
    dirs, direct_files = immediate_children(files, module_root)

    lines = [
        "## Repository explorer",
        "",
        START,
        "Generated from version-control candidates. Directories remain compact; "
        "use this module's curated sections for canonical files.",
        "",
    ]

    for filename in direct_files:
        lines.append(f"- {md_file_link(filename)}")
    if direct_files and dirs:
        lines.append("")
    for dirname in dirs:
        n = count_under(files, f"{module_root}/{dirname}")
        plural = "file" if n == 1 else "files"
        readme = f"{module_root}/{dirname}/README.md"
        if readme in files:
            lines.append(
                f"- **{esc(dirname)}/** — {n} {plural}; "
                f"{md_file_link(f'{dirname}/README.md', 'guide')}"
            )
        else:
            lines.append(f"- **{esc(dirname)}/** — {n} {plural}")
    lines += ["", END]

    return "\n".join(lines)

def replace_or_insert(text: str, generated: str) -> str:
    replacement = generated.rstrip() + "\n\n"

    if SECTION_RE.search(text):
        return SECTION_RE.sub(replacement, text, count=1)

    lines = text.splitlines()
    insert_at = None

    for index, line in enumerate(lines):
        if index > 0 and line.startswith("## "):
            insert_at = index
            break

    if insert_at is None:
        return text.rstrip() + "\n\n" + generated.rstrip() + "\n"

    before = "\n".join(lines[:insert_at]).rstrip()
    after = "\n".join(lines[insert_at:]).lstrip()

    return before + "\n\n" + generated.rstrip() + "\n\n" + after + "\n"

def update_file(path: Path, generated: str, check: bool) -> bool:
    if not path.exists():
        return False

    old = path.read_text(encoding="utf-8")
    new = replace_or_insert(old, generated)

    if new == old:
        return False

    if check:
        print(f"OUTDATED: {path}")
        return True

    path.write_text(new, encoding="utf-8")
    print(f"Updated: {path}")
    return True

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update collapsible repository navigation in README files."
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    root = repo_root()
    files = repository_files(root)
    changed = False

    changed |= update_file(
        root / "README.md",
        render_root(files),
        args.check,
    )

    for readme, module_root in MODULE_READMES.items():
        changed |= update_file(
            root / readme,
            render_module(files, module_root),
            args.check,
        )

    if args.check:
        if changed:
            print("Repository navigation is out of date.")
            return 1

        print("Repository navigation is up to date.")
        return 0

    if not changed:
        print("Repository navigation is already up to date.")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
