# Geodata documentation organization plan

This document defines how Geodata documentation is structured, how it maps to
`src/geodata`, and how we keep prose, notebooks, and API reference aligned as
the library evolves. It is intended for contributors working on the
**documentation branch** and for anyone opening a PR that changes user-facing
behavior.

**Status:** living plan (update this file when conventions change).

---

## 1. Goals

1. **One clear user journey** — readers should know whether to use the modern
   (`load_dataset` → model → optional masking) or legacy (`Dataset` → `Cutout` →
   `convert`) workflow without reading the entire site.
2. **Docs follow code** — every public API change in `src/` has a defined doc
   touchpoint (prose, notebook, or autoapi docstring).
3. **No orphan pages** — every `.md`, `.rst`, and `.ipynb` under
   `docs/source/` appears in a `toctree` or is explicitly marked as internal
   (see [Section 5](#5-file-types-and-conventions)).
4. **Reproducible examples** — tutorials should run offline where possible
   (ERA5 `*_test` fixtures) so CI and local builds do not depend on CDS
   credentials.
5. **Separation of concerns** — migration plans and design notes stay in
   `development/` or clearly labeled plan pages; user-facing tutorials stay
   task-focused.

---

## 2. Current state (baseline)

### 2.1 Two parallel workflows

Geodata currently exposes two stacks. Both are valid; documentation must label
them explicitly.

| Aspect | Modern workflow | Legacy workflow |
|--------|-----------------|-----------------|
| Data access | `geodata.datasets.load_dataset(...)` | `geodata.Dataset(module=..., weather_data_config=...)` |
| Subsetting | `BaseDataset` bounds / model `xs`/`ys` | `geodata.Cutout` + `prepare()` |
| Transform | `geodata.model.wind`, `geodata.model.pvlib` | `geodata.convert.*` on Cutouts |
| Masking (apply) | `geodata.XarrayMask` | `cutout.add_mask()` + `cutout.mask()` |
| Masking (create) | `geodata.Mask` (same for both) | `geodata.Mask` (same for both) |
| Primary docs | `datasets/`, `modeling/` | `intro.rst`, mask Cutout notebooks |

**Canonical path for new features:** modern workflow. Legacy paths remain
documented until explicitly deprecated.

### 2.2 Documentation build stack

| Piece | Location | Role |
|-------|----------|------|
| Sphinx config | `docs/source/conf.py` | MyST, notebooks, autoapi |
| Site root | `docs/source/index.rst` | Top-level toctrees |
| Landing narrative | `docs/source/intro.rst` | Overview (still legacy-heavy) |
| API reference | autoapi → `src/geodata` | Generated from docstrings |
| Notebooks | `myst_nb`, `nb_execution_mode = "off"` | Committed outputs; not executed on build |

### 2.3 Known gaps (as of this plan)

| Gap | Impact | Priority | Status |
|-----|--------|----------|--------|
| `intro.rst` teaches legacy Cutout/convert as the main story | New users miss models + `XarrayMask` | P0 | **Done** — modern intro on homepage; legacy moved to `legacy/workflow.rst` |
| Modeling pages missing recent API options (`compact_output`, flexible `xs`/`ys`) | Docs diverge from `src` | P0 | **Done** — see modeling/wind/index and modeling/pvlib/index |
| Wind capacity-factor internals not in wind toctree | Deep-dive exists only in source/comments | P1 | **Done (Option A)** — “Understanding the output” in interpolation/extrapolation Step 5 |
| `xarray_mask_tutorial.ipynb` referenced by `xarray_mask_workflow.rst` but may be missing from tree | Broken `:doc:` link | P0 | **Done** — notebook added under `mask/` |
| Mask section mixes user tutorials with `mask_xarray_migration_plan.md` | Hard to tell “how-to” vs “plan” | P1 | **Done** — plans moved to `development/` |
| `development/offline-era5-fixture-datasets.md` not linked from modeling tutorials | Readers assume CDS required | P1 | **Done** — overview + era5 + intro link fixtures |
| Example scripts in `docs/source/mask/*.py` not classified | Unclear if maintained or one-off | P2 | Open |
| README points to placeholder doc URL | External discoverability | P2 | Open |

---

## 3. Target information architecture

Organize the site by **user task**, not by file type. Recommended sidebar
structure (matches `index.rst` with clearer intent):

```
Geodata docs
├── Getting started
│   ├── Package setup
│   ├── Supported I/O formats
│   └── Workflow chooser (NEW — short page: modern vs legacy)
├── Datasets
│   ├── Overview (load_dataset, list_datasets)
│   ├── ERA5 (CDS setup + configs)
│   ├── MERRA2
│   └── Weather data config reference
├── Modeling
│   ├── Wind (index + interpolation + extrapolation; CF notes in Step 5)
│   └── PVLib (index + future subpages)
├── Masking
│   ├── Create masks (mask_creation_workflow.ipynb)
│   ├── Apply with XarrayMask (xarray_mask_tutorial.ipynb)
│   ├── Troubleshoot (mask_troubleshoot.md)
│   └── Cutout apply → legacy/mask_on_cutout
├── Visualization
├── Development (contributors)
│   ├── Documentation organization (this file)
│   ├── Offline ERA5 fixtures
│   ├── mask_xarray_migration_plan
│   └── xarray_mask_workflow (implementation notes)
└── API reference (autoapi)
```

### 3.1 Page roles (Diátaxis)

Use four doc types consistently:

| Type | Purpose | Examples |
|------|---------|----------|
| **Tutorial** | Learning-oriented, step-by-step | Notebooks, `modeling/wind/interpolation.rst` |
| **How-to guide** | Goal-oriented recipe | `xarray_mask_tutorial.ipynb`, ERA5 CDS setup |
| **Reference** | Accurate, complete | autoapi, `weather_data_config.md`, turbine YAML lists |
| **Explanation** | Concepts and design | migration plans, wind CF summary in interpolation Step 5 |

Label migration/plan documents at the top:

```markdown
> **Audience:** contributors and maintainers. For usage, see [XarrayMask tutorial](../mask/xarray_mask_tutorial.ipynb).
```

---

## 4. Source code ↔ documentation map

Maintain this table when adding modules. **Primary doc** is the page that must
be updated first when behavior changes.

| `src/geodata` area | Primary doc | Secondary / API |
|--------------------|-------------|-----------------|
| `datasets/_base.py`, `datasets/era5/*` | `datasets/overview.rst`, `datasets/era5.rst` | autoapi |
| `datasets/merra2/*` (legacy) | `legacy/merra2/*` | autoapi |
| `datasets/era5/fixture.py` (`*_test`) | `development/offline-era5-fixture-datasets.md` | modeling tutorials (offline note) |
| `model/wind/*` | `modeling/wind/index.rst`, `interpolation.rst`, `extrapolation.rst` | autoapi |
| `model/pvlib/_base.py` | `modeling/pvlib/index.rst` | autoapi |
| `model/_base.py` (slice sel, I/O) | modeling pages (bounding box sections) | autoapi |
| `mask.py` (legacy Mask) | `mask/mask_creation_workflow.ipynb` | autoapi |
| `mask/xarray_mask.py`, `mask/spatial.py` | `mask/xarray_mask_tutorial.ipynb` | `development/xarray_mask_workflow.rst` (contributors) |
| `cutout.py`, `convert.py`, `preparation.py` | `legacy/workflow.rst` | autoapi |
| Cutout-based masking | `legacy/mask_on_cutout.ipynb` | — |
| `plot.py` | `visualization/visualization.ipynb` | autoapi |
| `resource.py`, `resources/*` | modeling pages (turbine/panel names) | — |
| `config.py` | `quick_start/packagesetup.md` | — |

### 4.1 Public exports (`__init__.py`)

When adding or removing symbols from `geodata.__all__`:

1. Update docstrings (autoapi).
2. Update `intro.rst` or the relevant tutorial if the symbol is part of a
   documented workflow.
3. Add a line to the [changelog section](#72-changelog-expectations) of the PR.

---

## 5. File types and conventions

### 5.1 Where files live

| Path | Use for |
|------|---------|
| `docs/source/quick_start/` | Install, env vars, I/O formats |
| `docs/source/datasets/` | Download, configs, dataset-specific outputs |
| `docs/source/modeling/<domain>/` | Model tutorials and domain index |
| `docs/source/mask/` | Mask tutorials, workflows, troubleshooting |
| `docs/source/visualization/` | Plotting notebooks |
| `docs/source/development/` | Contributor docs, fixtures, **this plan**, internal design |
| `docs/source/_static/` | Images referenced from rst/md |

### 5.2 Format choice

| Format | When to use |
|--------|-------------|
| `.rst` | Sphinx-native pages with toctrees (section indexes) |
| `.md` (MyST) | Prose guides, plans, troubleshooting |
| `.ipynb` | Executable narratives with plots; keep outputs committed |

### 5.3 Naming

- User-facing: `snake_case` or `kebab-case` descriptive names
  (`xarray_mask_tutorial.ipynb`, `mask_troubleshoot.md`).
- Plans: suffix or folder under `development/` (`*_plan.md`, `*_known_issues.md`).
- Example scripts: `docs/source/<topic>/examples/` (proposed) — not mixed with
  built pages unless listed in toctree.

### 5.4 Internal vs published pages

Pages under `development/` and `mask/*_plan.md` are **contributor-facing**.
They stay in the toctree under **Development** or with an audience banner so
users are not sent to migration checklists by mistake.

Optional future convention: prefix internal-only files with `_` and exclude in
`conf.py` `exclude_patterns` — not required if audience banners are used.

---

## 6. Keeping documentation up to date

### 6.1 PR checklist (code changes)

Every PR that changes `src/geodata` should answer:

- [ ] Does this change **public API** or default behavior?
- [ ] Which **primary doc** row in [Section 4](#4-source-code--documentation-map) applies?
- [ ] Are **docstrings** updated for autoapi?
- [ ] Is there a **minimal code snippet** in prose docs or a test that can be copied?
- [ ] Do **notebooks** need re-run outputs (if affected)?
- [ ] Does `intro.rst` need a **workflow label** (modern vs legacy) if touched?

If the answer to the first question is yes and no doc file is updated, the PR
should either include doc updates or link a follow-up issue.

### 6.2 Documentation-only PRs (this branch)

Recommended batching for the documentation branch:

| Phase | Work | Outcome |
|-------|------|---------|
| **A — Structure** | Workflow chooser; relabel legacy in `intro.rst`; wire orphan pages into toctrees | Clear navigation |
| **B — Sync with recent `src`** | `XarrayMask`, `compact_output`, slice/bounds notes, fixture offline path | Factual parity with code |
| **C — Depth** | Wind CF explanation, mask troubleshooting, merge-layer known issues | Explanation layer |
| **D — Hygiene** | Move example `.py` to `examples/`; README doc URL; trim stale “planned” notes | Lower maintenance cost |

### 6.3 When to update which layer

| Change in `src` | Update prose/notebook | Update docstrings only |
|-----------------|----------------------|-------------------------|
| New public class or method | Yes | Yes |
| New optional parameter with non-obvious default | Yes (one example) | Yes |
| Internal refactor, same API | No | Only if signatures changed |
| Bug fix affecting coordinates, units, or outputs | Yes (note in troubleshooting or tutorial) | Yes |
| New `*_test` fixture config | `development/offline-era5-fixture-datasets.md` | — |
| Deprecation | Yes + migration plan | Yes |

### 6.4 Single source of truth

| Content | Source of truth | Docs should… |
|---------|-----------------|--------------|
| Function signatures | `src/` + autoapi | Not duplicate parameter lists |
| End-to-end workflows | Notebooks + rst tutorials | Link to tests under `tests/pr/` |
| Dataset registry names | `list_datasets()` / `datasets/registry` | Regenerate or manually sync lists in `overview.rst` when configs added |
| Turbine/panel names | `resources/windturbine/`, `resources/solarpanel/` | Show representative examples, not full catalogs |

Prefer **short examples from tests** over hand-written snippets that drift:

```python
# Pattern: tests/pr/test_xarray_mask.py → docs/source/mask/xarray_mask_tutorial.ipynb
```

### 6.5 Build and review

Local build:

```bash
cd docs && make html
# open _build/html/index.html
```

Before merging the documentation branch:

1. `make html` completes without warnings for missing `:doc:` references.
2. New pages appear in the correct toctree (sidebar).
3. Notebooks render (committed outputs present; `nb_execution_mode` is `off`).
4. autoapi pages generate for new modules.

Future CI enhancements (optional):

- Sphinx `-W` (warnings as errors) on PRs touching `docs/`.
- Link check for internal `:doc:` and relative md links.
- Script to diff `list_datasets()` output against `overview.rst` mentions.

---

## 7. Immediate backlog for the documentation branch

Actionable items in recommended order.

### P0 — Navigation and broken links

1. ~~**Add workflow chooser**~~ — **Done:** homepage (`intro.rst`) is the modern workflow; legacy content lives under **Legacy workflow** (`legacy/workflow.rst`).
2. ~~**Ensure `xarray_mask_tutorial.ipynb` exists**~~ — **Done** (`docs/source/mask/xarray_mask_tutorial.ipynb`, included via mask `*` toctree).
3. ~~**Update `intro.rst` masking section**~~ — **Done:** modern intro uses `XarrayMask`; Cutout masking unchanged in `legacy/workflow.rst`.

### P0 — Sync with recent source changes

4. ~~**`modeling/pvlib/index.rst`** — document `compact_output`~~ — **Done**.
5. ~~**`modeling/wind/index.rst` and interpolation.rst** — document flexible `xs`/`ys`~~ — **Done**.
6. ~~**`datasets/era5.rst`** — clarify CDS download vs offline fixtures~~ — **Done:** overview has recommended download; era5.rst covers CDS setup + optional cdsapi verify.

### P1 — Structure and depth

7. ~~**Wind CF documentation**~~ — **Done (Option A):** expanded Step 5 in interpolation/extrapolation; no separate internals page.
8. ~~**Reorganize mask toctree intent**~~ — **Done:** Mask = creation + apply tutorials + troubleshoot; plans under Development.
9. **`mask/merge_layer_known_issues.md`** — publish under mask with troubleshooting cross-links.
10. **Link fixture doc from modeling tutorials** — one paragraph + code using `load_dataset("wind_3d_hourly_test")`.

### P2 — Hygiene

11. Move `docs/source/mask/create_mask.py`, `split_china.py`, etc. to `docs/source/mask/examples/` (exclude from glob toctree or document as examples).
12. Fix README documentation URL placeholder.
13. Audit `input_output.md` “planned” notes against `list_datasets()`.
14. Add **this plan** to `index.rst` Development toctree (done when this file is merged).

---

## 7.2 Changelog expectations

Documentation PRs should summarize:

- **User-visible:** what readers can now do or what corrected behavior is documented.
- **Structural:** new pages, moved pages, deprecated paths.
- **Not required:** typo fixes only.

For paired code+doc releases, use a single changelog entry covering both.

---

## 8. Long-term governance

### 8.1 Ownership (suggested)

| Area | Default maintainer focus |
|------|--------------------------|
| Datasets / ERA5 fixtures | Whoever changes `datasets/era5/` |
| Wind / PV modeling | Model module authors |
| Mask / XarrayMask | Mask package authors |
| Legacy Cutout/convert | Touch only when behavior changes; avoid new features here |

### 8.2 Deprecation policy for docs

When deprecating APIs:

1. Mark in docstring + autoapi.
2. Add “Deprecated” admonition in legacy tutorial.
3. Record timeline in a `development/` plan or release notes.
4. Remove legacy tutorial sections only after code removal or major version bump.

### 8.3 Quarterly doc audit (lightweight)

Every ~3 months or before a release:

1. Run `list_datasets()` and compare to `datasets/overview.rst`.
2. Scan `intro.rst` for legacy-only examples without modern pointers.
3. Grep docs for `planned`, `TODO`, `FIXME`.
4. Confirm `make html` clean build.
5. Update [Section 2.3](#23-known-gaps-as-of-this-plan) gap table in this file.

### 8.4 Relationship to API reference

autoapi is the **reference layer**; tutorials should not duplicate every
argument. Convention:

- Tutorials: one worked example with common options.
- Reference: full signatures via docstrings (NumPy style, `sphinx.ext.napoleon`).
- Explanation pages: algorithms and data flow (e.g. wind CF pipeline).

When adding a feature, **docstring first**, then **one tutorial paragraph** —
not a third full copy in markdown.

---

## 9. Appendix: proposed `index.rst` Development toctree

```rst
.. toctree::
   :maxdepth: 1
   :caption: Development
   :hidden:

   development/documentation-organization-plan
   development/offline-era5-fixture-datasets
```

Mask migration plan remains under `mask/` glob but should use a contributor
banner (see [Section 3.1](#31-page-roles-diátaxis)).

---

## 10. Appendix: doc branch merge strategy

1. **Land structure first** (toctrees, workflow chooser, intro labels) so follow-up edits have a home.
2. **Land content sync** (modeling/mask/datasets factual updates) in the same branch or stacked PRs by area.
3. **Avoid** mixing large narrative rewrites with unrelated code changes — keeps review focused.
4. After merge, tag a docs release note listing: new XarrayMask path, fixture-based tutorials, deprecated/legacy labeling.

---

## Document history

| Date | Change |
|------|--------|
| 2026-06-02 | Initial organization plan for documentation branch |
