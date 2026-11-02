"""Helpers for executing the Python programs embedded in rMAP_TB.wdl outside Cromwell.

Each heavy task carries its logic in a Python here-document inside the WDL
`command` block. These helpers extract that program for a named task,
substitute the WDL `~{...}` placeholders with test values, and run it in a
scratch directory with the local interpreter, so the exact code that ships in
the workflow is what gets tested.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WDL_PATH = REPO_ROOT / "rMAP_TB.wdl"

_PLACEHOLDER = re.compile(r"~\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}")
_IF_DEFINED = re.compile(r"^if\s+defined\((\w+)\)\s+then\s+(?:select_first\(\[\s*(\w+)\s*\]\)|(\w+))\s+else\s+\"\"$")
_DEFAULT = re.compile(r"^default=\"([^\"]*)\"\s+(\w+)$")


def wdl_text() -> str:
    return WDL_PATH.read_text(encoding="utf-8")


def task_body(name: str, text: str | None = None) -> str:
    text = text or wdl_text()
    headers = list(re.finditer(r"^(task|workflow)\s+(\w+)\s*\{", text, re.M))
    for i, h in enumerate(headers):
        if h.group(1) == "task" and h.group(2) == name:
            end = headers[i + 1].start() if i + 1 < len(headers) else len(text)
            return text[h.end():end]
    raise KeyError(f"task {name} not found in {WDL_PATH}")


def command_block(name: str) -> str:
    m = re.search(r"command\s*<<<(.*?)\n\s*>>>", task_body(name), re.S)
    if not m:
        raise KeyError(f"task {name} has no command block")
    return m.group(1)


def python_heredoc(name: str, delimiter: str = "PY") -> str:
    """Return the body of the `<<'DELIM' ... DELIM` here-document of a task."""
    cmd = command_block(name)
    m = re.search(rf"<<-?'?{delimiter}'?\n(.*?)\n{delimiter}(?:\n|$)", cmd, re.S)
    if not m:
        raise KeyError(f"task {name} has no {delimiter} here-document")
    return m.group(1)


def substitute(src: str, values: dict) -> str:
    """Replace `~{...}` placeholders with values from `values`.

    Supports plain names (`~{width}`), `~{if defined(x) then x else ""}` and
    `~{default="" x}`. Missing names resolve to the empty string.
    """
    def repl(m):
        expr = m.group(1).strip()
        im = _IF_DEFINED.match(expr)
        if im:
            key = im.group(1)
            return str(values.get(key, "") or "")
        dm = _DEFAULT.match(expr)
        if dm:
            return str(values.get(dm.group(2), dm.group(1)))
        if re.fullmatch(r"\w+", expr):
            return str(values.get(expr, ""))
        raise ValueError(f"unsupported placeholder: ~{{{expr}}}")
    out = _PLACEHOLDER.sub(repl, src)
    assert "~{" not in out, "unsubstituted placeholder remains"
    return out


def run_python(src: str, cwd: Path, args: list[str] | None = None, timeout: int = 300) -> subprocess.CompletedProcess:
    script = cwd / "_wdl_block.py"
    script.write_text(src)
    return subprocess.run([sys.executable, "-I", str(script), *(args or [])],
                          cwd=cwd, capture_output=True, text=True, timeout=timeout)
