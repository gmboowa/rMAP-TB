#!/usr/bin/env python3
"""Static checks for rMAP_TB.wdl that do not need a WDL engine.

The workflow embeds large Bash and Python programs inside `command <<< >>>`
blocks. A syntax error there is only discovered when the task runs, often
hours into a workflow. This script catches those errors in seconds:

  1. every `command <<< ... >>>` body passes `bash -n` after WDL placeholders
     are substituted;
  2. every Python here-document (`<<'PY...'`) compiles with the local Python;
  3. f-strings in those blocks are compatible with Python 3.11, the version in
     the `python:3.11-slim` task image (no re-use of the enclosing quote
     character or backslashes inside `{}` expressions);
  4. every `call` input refers to an existing task output, and every workflow
     output refers to an existing task output;
  5. no task declares an input that the call block passes but the task lacks.

Exit status is non-zero when any check fails. Run with `--verbose` for a
per-task summary.

Example:
    python3 scripts/check_wdl.py rMAP_TB.wdl
"""

from __future__ import annotations

import argparse
import ast
import io
import os
import re
import subprocess
import sys
import tempfile
import tokenize
from pathlib import Path

PLACEHOLDER = re.compile(r"~\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}")
PY_HEREDOC = re.compile(r"<<-?'?(PY[A-Z_]*)'?\n(.*?)\n\1(?:\n|$)", re.S)


def substitute(body: str) -> str:
    return PLACEHOLDER.sub("X", body)


def fstring_py311_issues(src: str) -> list[tuple[int, int]]:
    """Positions of f-string constructs that Python 3.11 rejects."""
    problems = []
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(src).readline))
    except Exception:
        return problems
    if not hasattr(tokenize, "FSTRING_START"):
        return problems  # local Python < 3.12 would already have failed to compile
    stack = []
    for t in toks:
        if t.type == tokenize.FSTRING_START:
            stack.append(t.string.lstrip("rRfFbB"))
        elif t.type == tokenize.FSTRING_END:
            stack.pop()
        elif stack and t.type != tokenize.FSTRING_MIDDLE:
            q = stack[-1]
            if t.type == tokenize.STRING and (q in t.string or "\\" in t.string):
                problems.append(t.start)
            elif t.type == tokenize.FSTRING_START and q in t.string.lstrip("rRfFbB"):
                problems.append(t.start)
            elif t.type not in (tokenize.STRING, tokenize.FSTRING_START) and q in t.string:
                problems.append(t.start)
    return problems


def split_tasks(text: str) -> dict[str, str]:
    """Split on top-level `task`/`workflow` headers. A body runs to the next
    header; this is robust to `}` at column 0 inside embedded Python."""
    headers = list(re.finditer(r"^(task|workflow)\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{", text, re.M))
    tasks = {}
    for i, h in enumerate(headers):
        if h.group(1) != "task":
            continue
        end = headers[i + 1].start() if i + 1 < len(headers) else len(text)
        tasks[h.group(2)] = text[h.end():end]
    return tasks


def after_command(body: str) -> str:
    m = re.search(r"command\s*<<<.*?\n\s*>>>", body, re.S)
    return body[m.end():] if m else body


def task_outputs(body: str) -> set[str]:
    m = re.search(r"\n\s*output\s*\{(.*?)\n\s*\}", after_command(body), re.S)
    if not m:
        return set()
    return set(re.findall(r"^\s*(?:Array\[[^\]]+\]|File|String|Int|Float|Boolean)\??\+?\s+([A-Za-z_][A-Za-z0-9_]*)\s*=", m.group(1), re.M))


def task_inputs(body: str) -> set[str]:
    m = re.search(r"\n\s*input\s*\{(.*?)\n\s*\}", body.split("command", 1)[0], re.S)
    if not m:
        return set()
    return set(re.findall(r"^\s*(?:Array\[[^\]]+\]|File|String|Int|Float|Boolean)\??\+?\s+([A-Za-z_][A-Za-z0-9_]*)", m.group(1), re.M))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("wdl", nargs="?", type=Path, default=Path(__file__).resolve().parents[1] / "rMAP_TB.wdl")
    ap.add_argument("--verbose", "-v", action="store_true")
    a = ap.parse_args(argv)

    text = a.wdl.read_text(encoding="utf-8")
    failures = 0
    tasks = split_tasks(text)
    if not tasks:
        print("ERROR: no tasks found")
        return 1

    for name, body in tasks.items():
        cm = re.search(r"command\s*<<<(.*?)\n\s*>>>", body, re.S)
        if not cm:
            print(f"ERROR  {name}: no command <<< >>> block")
            failures += 1
            continue
        cmd = substitute(cm.group(1))

        with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as fh:
            fh.write(cmd)
            tmp = fh.name
        try:
            r = subprocess.run(["bash", "-n", tmp], capture_output=True, text=True)
        finally:
            os.unlink(tmp)
        if r.returncode:
            failures += 1
            print(f"ERROR  {name}: bash syntax\n{r.stderr.strip()}")

        blocks = PY_HEREDOC.findall(cmd)
        py_issues = 0
        for delim, src in blocks:
            try:
                ast.parse(src)
            except SyntaxError as exc:
                failures += 1
                py_issues += 1
                print(f"ERROR  {name}: Python block {delim} does not compile: {exc}")
                continue
            bad = fstring_py311_issues(src)
            if bad:
                failures += 1
                py_issues += 1
                print(f"ERROR  {name}: Python block {delim} has {len(bad)} f-string construct(s) rejected by Python 3.11 at lines {sorted({l for l, _ in bad})[:10]}")
        if a.verbose:
            print(f"ok     {name}: bash ok, {len(blocks)} Python block(s), {py_issues} issue(s)")

    # Wiring checks on the workflow block
    wm = re.search(r"^workflow\s+\w+\s*\{", text, re.M)
    if not wm:
        print("ERROR: workflow block not found")
        return 1
    nxt = re.search(r"^task\s+\w+\s*\{", text[wm.end():], re.M)
    wf = text[wm.end(): wm.end() + nxt.start()] if nxt else text[wm.end():]
    outputs_by_task = {t: task_outputs(b) for t, b in tasks.items()}
    inputs_by_task = {t: task_inputs(b) for t, b in tasks.items()}

    for cm in re.finditer(r"call\s+(\w+)\s*\{\s*input:(.*?)\n\s*\}", wf, re.S):
        callee, block = cm.group(1), cm.group(2)
        if callee not in tasks:
            failures += 1
            print(f"ERROR  call {callee}: task not defined")
            continue
        for im in re.finditer(r"^\s*(\w+)\s*=\s*(.+?),?\s*$", block, re.M):
            iname, expr = im.group(1), im.group(2)
            if iname not in inputs_by_task[callee]:
                failures += 1
                print(f"ERROR  call {callee}: passes input {iname!r} that the task does not declare")
            for ref in re.findall(r"\b([A-Z][A-Z0-9_]+)\.([a-z_][a-z0-9_]*)\b", expr):
                t, o = ref
                if t in tasks and o not in outputs_by_task[t]:
                    failures += 1
                    print(f"ERROR  call {callee}: input {iname} references {t}.{o}, which {t} does not output")

    om = None
    for om in re.finditer(r"\n\s*output\s*\{(.*?)\n\s*\}", wf, re.S):
        pass
    if om:
        for t, o in re.findall(r"=\s*([A-Z][A-Z0-9_]+)\.([a-z_][a-z0-9_]*)", om.group(1)):
            if t in tasks and o not in outputs_by_task[t]:
                failures += 1
                print(f"ERROR  workflow output references {t}.{o}, which {t} does not output")

    print(f"{len(tasks)} task(s) checked; {failures} problem(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
