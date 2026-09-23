#!/usr/bin/env python3
"""Update deterministic repository trees in README files.

The tree is built from Git-tracked files (the index), not arbitrary local files.
This keeps generated outputs, temporary files and ignored artifacts out of public
README documentation.

Usage:
    python3 scripts/update_repository_tree.py
    python3 scripts/update_repository_tree.py --check
    python3 scripts/update_repository_tree.py --stage
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

START_MARKER = "<!-- AUTO_REPOSITORY_TREE_START -->"
END_MARKER = "<!-- AUTO_REPOSITORY_TREE_END -->"
SECTION_HEADING = "## Repository structure"


@dataclass(frozen=True)
class ReadmeTarget:
    path: str
    prefix: str = ""
    max_depth: int | None = None


TARGETS = (
    # Root: compact navigation. Detailed trees live in module READMEs.
    ReadmeTarget("README.md", "", 4),
    ReadmeTarget("alignment_analysis/README.md", "alignment_analysis/", None),
    ReadmeTarget("assemblers/README.md", "assemblers/", 4),
    ReadmeTarget("assemblers/whole_genome_asm/README.md", "assemblers/whole_genome_asm/", None),
    ReadmeTarget("assemblers/whole_genome_asm/assessment/README.md", "assemblers/whole_genome_asm/assessment/", None),
    ReadmeTarget("assembly_analysis/README.md", "assembly_analysis/", None),
    ReadmeTarget("variant_calling_analysis/README.md", "variant_calling_analysis/", None),
)


class Node(dict):
    """Simple nested-dict tree node."""


def run_git(repo_root: Path, *args: str) -> str:
    proc = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "git command failed")
    return proc.stdout


def find_repo_root() -> Path:
    proc = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if proc.returncode != 0:
        raise SystemExit("ERROR: run this command from inside the Git repository")
    return Path(proc.stdout.strip()).resolve()


def tracked_files(repo_root: Path) -> list[str]:
    # --cached reads the Git index. Newly staged files are therefore included
    # before the commit is created, which is exactly what the pre-commit hook
    # needs.
    output = run_git(repo_root, "ls-files", "--cached")
    return sorted(line.strip() for line in output.splitlines() if line.strip())


def relative_files(files: Iterable[str], prefix: str) -> list[str]:
    if not prefix:
        return list(files)
    return [path[len(prefix):] for path in files if path.startswith(prefix) and path != prefix.rstrip("/")]


def build_tree(paths: Iterable[str], max_depth: int | None) -> Node:
    root: Node = Node()
    truncated_dirs: set[tuple[str, ...]] = set()

    for raw_path in paths:
        parts = [part for part in Path(raw_path).parts if part not in {"", "."}]
        if not parts:
            continue

        if max_depth is not None and len(parts) > max_depth:
            truncated_dirs.add(tuple(parts[:max_depth]))
            parts = parts[:max_depth]

        node = root
        for part in parts:
            node = node.setdefault(part, Node())

    for parts in sorted(truncated_dirs):
        node = root
        for part in parts:
            node = node.setdefault(part, Node())
        node.setdefault("…", Node())

    return root


def sort_items(node: Node):
    # Directories first, files second; alphabetical within each class.
    def key(item):
        name, children = item
        is_dir = bool(children) or name == "…"
        return (0 if is_dir else 1, name.lower())

    return sorted(node.items(), key=key)


def render_tree(node: Node, root_label: str) -> str:
    lines = [root_label.rstrip("/") + "/"]

    def walk(current: Node, prefix: str) -> None:
        items = sort_items(current)
        for index, (name, children) in enumerate(items):
            last = index == len(items) - 1
            connector = "└── " if last else "├── "
            suffix = "/" if children and name != "…" else ""
            lines.append(prefix + connector + name + suffix)
            if children and name != "…":
                extension = "    " if last else "│   "
                walk(children, prefix + extension)

    walk(node, "")
    return "\n".join(lines)


def generated_block(tree_text: str, detail_note: str | None = None) -> str:
    note = (
        "This section is generated from Git-tracked files. "
        "Do not edit the tree manually."
    )
    if detail_note:
        note += " " + detail_note

    return (
        f"{START_MARKER}\n"
        f"{note}\n\n"
        f"```text\n{tree_text}\n```\n"
        f"{END_MARKER}"
    )


def insert_or_replace(readme_text: str, block: str) -> str:
    if START_MARKER in readme_text or END_MARKER in readme_text:
        if START_MARKER not in readme_text or END_MARKER not in readme_text:
            raise ValueError("README has only one auto-tree marker; repair markers before rerunning")
        start = readme_text.index(START_MARKER)
        end = readme_text.index(END_MARKER, start) + len(END_MARKER)
        return readme_text[:start] + block + readme_text[end:]

    section = f"{SECTION_HEADING}\n\n{block}\n\n"

    # Place the tree near the beginning: before the first H2 section. This
    # preserves the title/introduction while making repository navigation easy.
    lines = readme_text.splitlines(keepends=True)
    insert_at = None
    offset = 0
    for line in lines:
        if line.startswith("## "):
            insert_at = offset
            break
        offset += len(line)

    if insert_at is None:
        base = readme_text.rstrip()
        return base + "\n\n" + section if base else section

    before = readme_text[:insert_at].rstrip()
    after = readme_text[insert_at:].lstrip("\n")
    return before + "\n\n" + section + after


def update_target(repo_root: Path, files: list[str], target: ReadmeTarget, *, write: bool = True) -> tuple[Path, bool]:
    readme = repo_root / target.path
    if not readme.is_file():
        return readme, False

    rel = relative_files(files, target.prefix)
    tree = build_tree(rel, target.max_depth)
    root_label = Path(target.prefix.rstrip("/")).name if target.prefix else repo_root.name

    detail_note = None
    if target.path == "README.md":
        detail_note = "The root tree is compact; module READMEs contain more detailed trees."

    block = generated_block(render_tree(tree, root_label), detail_note)
    old = readme.read_text(encoding="utf-8")
    new = insert_or_replace(old, block)
    changed = new != old
    if changed and write:
        readme.write_text(new, encoding="utf-8")
    return readme, changed


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="fail if README trees are out of date")
    mode.add_argument("--stage", action="store_true", help="update trees and stage changed README files")
    args = parser.parse_args()

    repo_root = find_repo_root()
    files = tracked_files(repo_root)
    changed: list[Path] = []

    for target in TARGETS:
        readme, was_changed = update_target(repo_root, files, target, write=not args.check)
        if was_changed:
            changed.append(readme)

    if args.check:
        if changed:
            print("Repository tree documentation is out of date:")
            for path in changed:
                print("  " + str(path.relative_to(repo_root)))
            print("Run: python3 scripts/update_repository_tree.py")
            return 1
        print("Repository tree documentation is up to date.")
        return 0

    if changed:
        print("Updated repository tree documentation:")
        for path in changed:
            rel = str(path.relative_to(repo_root))
            print("  " + rel)
            if args.stage:
                run_git(repo_root, "add", "--", rel)
        if args.stage:
            print("Updated README files were staged.")
    else:
        print("Repository tree documentation is already up to date.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
