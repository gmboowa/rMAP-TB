# Tests

The suite exercises the programs embedded in `rMAP_TB.wdl` without Cromwell, Docker
or network access. `tests/wdl_blocks.py` extracts the Python here-document of a task
from the WDL, substitutes its `~{...}` placeholders, and runs it with the local
interpreter, so the code under test is exactly the code that ships in the workflow.

| File | Scope |
|---|---|
| `test_wdl_static.py` | Task inventory, `scripts/check_wdl.py` (Bash and Python syntax, Python 3.11 f-string compatibility, call wiring), regressions that must not return (1×1 placeholder PNG, csv field limit, heatmap title) |
| `test_tree_visualization.py` | `TREE_VISUALIZATION`: SVG tree with annotations, >131 KB TB-Profiler fields, placeholder raster legibility, Newick edge cases (quoted labels, star tree, zero lengths, reference pruning, 300 tips), IQ-TREE fallback handling, bootstrap labelling threshold |
| `test_surveillance_visuals.py` | `TB_SURVEILLANCE_SUMMARY_VISUALS`: heatmap without title and without clipped labels, log-scale legend, empty-matrix message, metadata and QC-rationale tables |
| `test_merge_report.py` | `MERGE_TB_REPORTS`: rendered tree embedded as SVG, renderer failure shown with log, legacy 68-byte PNG rejected, incomplete annotations flagged, coverage-gate and IQ-TREE exclusions explained, report self-contained |
| `test_scripts.py` | `make_inputs_json.py`, `validate_inputs.py`, `collect_outputs.py`, `check_wdl.py`, example input files |

Run from the repository root:

```bash
python3 -m pip install pytest
python3 -m pytest -q tests/
```

Requirements: Python 3.11 or later and Bash. ETE3 is optional; when it is absent the
tests verify the placeholder path of the raster renderer instead. The synthetic panel
used by the tests (13 MTBC isolates, 3 NTM, one identical pair, one distant outlier,
one mixed-infection flag) is generated in `conftest.py`.
