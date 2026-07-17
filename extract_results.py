#!/usr/bin/env python3

# This file is part of g-mm-repair
# <https://github.com/ftosoni/g-mm-repair>.
# Copyright (c) 2026 Francesco Tosoni.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Single extraction layer: raw experiment logs -> manuscript tables + figure data.

This is the *only* place where measured numbers cross from the logs into the paper.
It parses the canonical raw logs written by ``reproduce.sh`` (see REPRODUCIBILITY.md)
and emits:

  * ``manuscript/tables/tab_*.tex``      -- the tabular bodies ``\\input{}``-ed by main.tex
  * ``manuscript/figures/data/*.dat``    -- pgfplots data read by the standalone TikZ figures

Every table/figure has one generator function whose docstring names the log(s) it reads,
so the log->artifact mapping is auditable. A missing input log is skipped with a warning
(not an error), so the pipeline can be run incrementally as experiments land.

Canonical logs all use ``## <key>`` block separators (one per dataset / configuration);
see ``read_blocks``. Run with ``python3 extract_results.py`` from the repo root.
"""
import os
import re
import sys

ROOT      = os.path.dirname(os.path.abspath(__file__))
LOG_DIR   = os.path.join(ROOT, "manuscript", "logs")
TAB_DIR   = os.path.join(ROOT, "manuscript", "tables")
FIGDAT_DIR= os.path.join(ROOT, "manuscript", "figures", "data")

# --------------------------------------------------------------------------------------
# Dataset configuration. The single source of truth for display names, dimensions, and
# table ordering. `key` is the block marker used in the logs (## <key>); `tex` is the
# LaTeX row label reproduced verbatim from the manuscript so tables render identically.
# --------------------------------------------------------------------------------------
GENO_REAL = [
    ("geno22",     r"Chr22 ($2504\times0.10\text{M}$)", 2504, 100000,  12.6e6),
    ("geno22full", r"Chr22 ($2504\times1.06\text{M}$)", 2504, 1055454, 128.0e6),
    ("geno21",     r"Chr21 ($2504\times0.10\text{M}$)", 2504, 100000,  14.0e6),
    ("geno21full", r"Chr21 ($2504\times1.05\text{M}$)", 2504, 1054447, 140.2e6),
    ("geno20",     r"Chr20 ($2504\times0.10\text{M}$)", 2504, 100000,  13.5e6),
    ("geno20full", r"Chr20 ($2504\times1.74\text{M}$)", 2504, 1739315, 201.0e6),
]
GENO_SYNTH = [
    ("geno_synth_small",     r"$\text{synth\_small } (2\text{K}\times50\text{K})$",    2000,  50000,  16616951),
    ("geno_synth_large",     r"$\text{synth\_large } (5\text{K}\times200\text{K})$",   5000,  200000, 152440925),
    ("geno_synth_ld_high",   r"$\text{synth\_ld\_high } (5\text{K}\times100\text{K})$",5000,  100000, 73471493),
    ("geno_synth_ld_low",    r"$\text{synth\_ld\_low } (5\text{K}\times100\text{K})$", 5000,  100000, 77009034),
    ("geno_synth_ind_large", r"$\text{synth\_ind\_large } (10\text{K}\times50\text{K})$",10000,50000,  71602677),
]
CROSSOVER = ("crossover_synth", r"$\text{crossover\_synth } (10\text{K}\times700\text{K})$", 10000, 700000, 1002.97e6)
GENO = GENO_REAL + GENO_SYNTH

# Short row labels for the single-column structural table (tab:geno_through only);
# the wide tables (tab_geno_time/tab_geno/tab_build) keep the full labels with dims.
SHORT_GENO = {
    "geno22": "Chr22", "geno22full": "Chr22 full",
    "geno21": "Chr21", "geno21full": "Chr21 full",
    "geno20": "Chr20", "geno20full": "Chr20 full",
    "geno_synth_small":     r"\texttt{synth\_small}",
    "geno_synth_large":     r"\texttt{synth\_large}",
    "geno_synth_ld_high":   r"\texttt{synth\_ld\_high}",
    "geno_synth_ld_low":    r"\texttt{synth\_ld\_low}",
    "geno_synth_ind_large": r"\texttt{synth\_ind\_large}",
}

# Wikidata relations (tab:graph_struct / tab:graph). tex is the \texttt row label.
# nnz (edge count) is a fixed dataset property, carried here (like GRAPH_SCALE) because
# the struct log does not emit it; it is also documented in REPRODUCIBILITY.md (§2.D).
GRAPH = [
    ("wd_sports_team", r"\texttt{wd\_sports\_team}", 332121,  29854,   1136249),
    ("wd_cast_member", r"\texttt{wd\_cast\_member}", 173977,  144095,  1033124),
    ("wd_citizenship", r"\texttt{wd\_citizenship}", 2874250, 2556,     3063058),
    ("wd_occupation",  r"\texttt{wd\_occupation}",  3459933, 10610,    4596658),
    ("wd_subclass_of", r"\texttt{wd\_subclass\_of}",1487709, 73417,    2024347),
]
GRAPH_SCALE = [
    ("wd_country",    r"\texttt{wd\_country} ($10.1$M)",    10058956, 1747,     10089284),
    ("wd_cites_work", r"\texttt{wd\_cites\_work} ($166.7$M)",7072574, 12245945, 166682725),
    ("swh_full",      r"\texttt{swh} ($1.22$G)",            45691499, 45691499, 1218488928),
]

REAL_HDR  = r"\multicolumn{%d}{l}{\textbf{Real Genotypes (1000 Genomes)}} \\"
SYNTH_HDR = r"\multicolumn{%d}{l}{\textbf{Synthetic Genotypes (Haplotypes)}} \\"

# --------------------------------------------------------------------------------------
# Generic log helpers
# --------------------------------------------------------------------------------------
def read_blocks(name):
    """Return {key: block_text} for a canonical ``## <key>`` separated log, or None.

    Everything before the first ``## `` marker is ignored (provenance header). Accepts
    a bare key on the marker line and stops the key at the first whitespace, so both
    ``## geno22`` and ``## geno22 (2504 x 100000)`` map to key ``geno22``.
    """
    path = os.path.join(LOG_DIR, name)
    if not os.path.exists(path):
        print(f"  [skip] missing log: manuscript/logs/{name}", file=sys.stderr)
        return None
    text = open(path, encoding="utf-8", errors="replace").read()
    blocks, key, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^##\s+(\S+)", line)
        if m:
            if key is not None:
                blocks[key] = "\n".join(buf)
            key, buf = m.group(1), []
        elif key is not None:
            buf.append(line)
    if key is not None:
        blocks[key] = "\n".join(buf)
    return blocks


def grabf(rx, text, default=None):
    m = re.search(rx, text)
    return float(m.group(1)) if m else default


def grabi(rx, text, default=None):
    m = re.search(rx, text)
    return int(m.group(1)) if m else default


def write_table(name, header_lines, rows, colspec):
    """Write a ``tabular`` body (no caption/label -- those stay in main.tex)."""
    os.makedirs(TAB_DIR, exist_ok=True)
    out = [r"% GENERATED by extract_results.py -- do not edit by hand.",
           r"\begin{tabular}{%s}" % colspec, r"\toprule"]
    out += header_lines
    out.append(r"\midrule")
    out += rows
    out += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB_DIR, name), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"  [ok]  manuscript/tables/{name}")


def write_dat(name, header, rows):
    os.makedirs(FIGDAT_DIR, exist_ok=True)
    out = ["# GENERATED by extract_results.py -- do not edit by hand.", header]
    out += rows
    with open(os.path.join(FIGDAT_DIR, name), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"  [ok]  manuscript/figures/data/{name}")


# number formatting -------------------------------------------------------------
def abbr(v):
    """Genotype-style magnitude abbreviation: M with 2 decimals, K with 1."""
    if v is None:
        return "--"
    v = float(v)
    if v >= 1e6:
        return f"{v/1e6:.2f}M"
    if v >= 1e3:
        return f"{v/1e3:.1f}K"
    return f"{v:.0f}"


def grp(v):
    """Comma-grouped integer (graph-table style)."""
    return "--" if v is None else f"{int(round(v)):,}"


def f2(v):  return "--" if v is None else f"{v:.2f}"
def f3(v):  return "--" if v is None else f"{v:.3f}"
def mb(b):  return "--" if b is None else f"{b/1e6:.1f}"   # bytes -> MB, 1 decimal


# per-block structural parse (shared by struct tables) --------------------------
STRUCT_RX = r"Max depth \(L\): (\d+), Max width \(w\*\): (\d+), Total layered rules: (\d+) \(raw NTs: (\d+), \+pt: (\d+)\)"

def parse_struct(block):
    m = re.search(STRUCT_RX, block)
    if not m:
        return None
    L, w, total, raw, pt = map(int, m.groups())
    return dict(L=L, w=w, total=total, raw=raw, pt=pt)


# --------------------------------------------------------------------------------------
# Table 1 -- tab:geno_through  (genotype structural figures)
#   log: logs/geno_struct.log   (gpu_test structural output, one ## <key> block each)
# --------------------------------------------------------------------------------------
def gen_tab_geno_through():
    # Structural lines also appear in the space/energy log; use it as a fallback.
    blocks = read_blocks("geno_struct.log") or read_blocks("geno_space_energy.log")
    if not blocks:
        return
    rows = [REAL_HDR % 5]
    for group, hdr in ((GENO_REAL, None), (GENO_SYNTH, SYNTH_HDR % 5)):
        if hdr:
            rows.append(hdr)
        for key, tex, r, c, nnz in group:
            s = parse_struct(blocks.get(key, ""))
            if not s:
                print(f"  [warn] no struct for {key}", file=sys.stderr)
                continue
            label = SHORT_GENO.get(key, tex)   # short labels: single-column table
            rows.append(f"{label} & {abbr(s['raw'])} & {s['L']} & {abbr(s['w'])} & {abbr(s['pt'])}\\\\")
    write_table("tab_geno_through.tex",
                [r"matrix & $|\mathcal{R}|$ & $L$ & $w^{*}$ & $+\text{pt}$\\"],
                rows, "lcccc")


# --------------------------------------------------------------------------------------
# Table 2 -- tab:geno_time  (average ms/vector)
#   logs: manuscript/logs/bio_results.csv     (GPU,OMP,CPU_Seq,mm_seq,mm_16)
#         manuscript/logs/geno_space_energy.log  (cuSPARSE ms; see gen_tab_geno)
# --------------------------------------------------------------------------------------
def _csv_map():
    path = os.path.join(LOG_DIR, "bio_results.csv")
    if not os.path.exists(path):
        print("  [skip] missing log: manuscript/logs/bio_results.csv", file=sys.stderr)
        return None
    out = {}
    for i, line in enumerate(open(path, encoding="utf-8")):
        if i == 0 or not line.strip():
            continue
        p = line.strip().split(",")
        num = lambda x: float(x) if x not in ("", "N/A") else None
        out[p[0]] = dict(gpu=num(p[1]), omp=num(p[2]), seq=num(p[3]),
                         mm_seq=num(p[4]), mm_16=num(p[5]))
    return out


def gen_tab_geno_time():
    csv = _csv_map()
    se  = _space_energy_map()   # for the cuSPARSE ms column
    if not csv:
        return
    rows = []
    for group, hdr in ((GENO_REAL, REAL_HDR % 8), (GENO_SYNTH, SYNTH_HDR % 8)):
        rows.append(hdr)
        for key, tex, r, c, nnz in group:
            d = csv.get(key)
            if not d:
                print(f"  [warn] no time row for {key}", file=sys.stderr)
                continue
            cus = (se or {}).get(key, {}).get("cus_ms")
            xseq = (d["mm_seq"] / d["gpu"]) if (d["mm_seq"] and d["gpu"]) else None
            rows.append(
                f"{tex} & {f2(d['gpu'])} & {f2(d['omp'])} & {f2(d['seq'])} & "
                f"{f2(d['mm_seq'])} & {f2(d['mm_16'])} & {f2(cus)} & "
                f"{('%.1f' % xseq) if xseq else '--'}\\\\")
    write_table("tab_geno_time.tex",
                [r"matrix & GPU & OpenMP & seq. & mmr (seq) & mmr (16th) & cuSPARSE & $\times$seq\\"],
                rows, "lrrrrrrr")


# --------------------------------------------------------------------------------------
# Table 3 -- tab:geno  (space & energy vs cuSPARSE, incl. the crossover_synth scale point)
# Figure 3 -- fig:space (footprint vs nnz)
#   log: manuscript/logs/geno_space_energy.log
#        one ## <key> block per dataset, containing gpu_test raw then cusparse_test raw.
#        Engine lines:  "GPU kernel: ... = X ms/vector", "MEM analytic bytes: N",
#                       "MEM peak bytes: N", "ENERGY mJ/vector: X"
#        cuSPARSE lines: "Average time: X ms/vector", "MEM analytic bytes: N" (or OOM)
# --------------------------------------------------------------------------------------
def _space_energy_map():
    blocks = read_blocks("geno_space_energy.log")
    if not blocks:
        return None
    out = {}
    for key, b in blocks.items():
        # Split engine vs cuSPARSE halves. reproduce.sh emits an explicit
        # "--- CUSPARSE ---" marker; fall back to a heuristic banner match.
        if "--- CUSPARSE ---" in b:
            eng_txt, cus_txt = b.split("--- CUSPARSE ---", 1)
        else:
            parts = re.split(r"(?im)^.*cusparse.*(?:baseline|test|CSR).*$", b, maxsplit=1)
            eng_txt = parts[0]
            cus_txt = parts[1] if len(parts) > 1 else b
        oom = bool(re.search(r"out of memory|OOM", cus_txt, re.I))
        out[key] = dict(
            eng_ms   = grabf(r"GPU kernel:.*=\s*([\d.]+) ms/vector", eng_txt),
            eng_mb   = (grabi(r"MEM analytic bytes:\s*(\d+)", eng_txt) or 0) / 1e6 or None,
            eng_pkmb = (grabi(r"MEM peak bytes:\s*(\d+)", eng_txt) or 0) / 1e6 or None,
            eng_mj   = grabf(r"ENERGY mJ/vector:\s*([\d.eE+-]+)", eng_txt),
            # nnz as counted by the reconstructed-CSR builder in the log; the per-dataset
            # constant above is only a fallback for blocks with no cuSPARSE half.
            nnz      = grabi(r"NNZ:\s*(\d+)", cus_txt),
            cus_ms   = None if oom else grabf(r"Average time:\s*([\d.]+) ms/vector", cus_txt),
            cus_mb   = None if oom else ((grabi(r"MEM analytic bytes:\s*(\d+)", cus_txt) or 0)/1e6 or None),
            cus_mj   = None if oom else grabf(r"ENERGY mJ/vector:\s*([\d.eE+-]+)", cus_txt),
            oom      = oom,
        )
    return out


def gen_tab_geno():
    se = _space_energy_map()
    if not se:
        return
    rows = []
    for group, hdr in ((GENO_REAL, REAL_HDR % 8), (GENO_SYNTH, SYNTH_HDR % 8)):
        rows.append(hdr)
        for key, tex, r, c, nnz in group:
            d = se.get(key)
            if not d:
                print(f"  [warn] no space/energy for {key}", file=sys.stderr)
                continue
            rows.append(
                f"{tex} & {abbr(d['nnz'] or nnz)} & {mb(d['eng_mb']*1e6 if d['eng_mb'] else None)} & "
                f"{mb(d['cus_mb']*1e6 if d['cus_mb'] else None)} & {f2(d['eng_ms'])} & "
                f"{f2(d['cus_ms'])} & {_i(d['eng_mj'])} & {_i(d['cus_mj'])}\\\\")
    # crossover row: at this scale cuSPARSE's CSR now fits and runs; only fall back to
    # "OOM" if the log actually records an out-of-memory cuSPARSE run.
    key, tex, r, c, nnz = CROSSOVER
    d = se.get(key)
    if d:
        cus_mb = "OOM" if d['oom'] else mb(d['cus_mb']*1e6 if d['cus_mb'] else None)
        cus_ms = "OOM" if d['oom'] else f2(d['cus_ms'])
        cus_mj = "OOM" if d['oom'] else _i(d['cus_mj'])
        rows.append(
            f"{tex} & {abbr(d['nnz'] or nnz)} & {mb(d['eng_mb']*1e6 if d['eng_mb'] else None)} & {cus_mb} & "
            f"{f2(d['eng_ms'])} & {cus_ms} & {_i(d['eng_mj'])} & {cus_mj}\\\\")
    header = [
        r"& & \multicolumn{2}{c}{space (MB)} & \multicolumn{2}{c}{time (ms)} & \multicolumn{2}{c}{energy (mJ)}\\",
        r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}",
        r"matrix ($\text{rows}\times\text{cols}$) & nnz & eng. & cuS. & eng. & cuS. & eng. & cuS.\\",
    ]
    write_table("tab_geno.tex", header, rows, "lr rr rr rr")


def _i(v):  # energy shown as integer mJ in the manuscript
    return "--" if v is None else f"{v:.0f}"


# --------------------------------------------------------------------------------------
# Table B -- tab:build  (two-part construction cost -- sec:limitations (ii))
# Two distinct, additive host costs are reported per matrix:
#   * grammar (RePair) build   -- the offline mm-repair compressor cost, seconds, SHARED
#       with the CPU baseline. log: manuscript/logs/grammar_build.log, matrepair "Compression Report"
#       (">> total time: X"), one ## <key> block each.
#   * +pt schedule build       -- the engine's MARGINAL host step (pass-through completion
#       + level-bucketing + terminal compaction), milliseconds. log: geno_space_energy.log,
#       "BUILD schedule ms: X" (or legacy "Schedule builder construction time: X ms").
# The amortization threshold (>= how many right-products before the build pays off) is
# computed on the marginal +pt cost only: build_ms/eval_ms. The grammar cost is the
# one-time price of the compressed format, amortized over the matrix's stored lifetime.
# --------------------------------------------------------------------------------------
def gen_tab_build():
    se  = read_blocks("geno_space_energy.log")
    gr  = read_blocks("grammar_build.log") or {}
    if not se:
        return
    rows = []
    for group, hdr in ((GENO_REAL, REAL_HDR % 4), (GENO_SYNTH, SYNTH_HDR % 4)):
        rows.append(hdr)
        for key, tex, r, c, nnz in group:
            eng = se.get(key, "").split("--- CUSPARSE ---", 1)[0]
            build_ms = grabf(r"BUILD schedule ms:\s*([\d.]+)", eng) \
                       or grabf(r"Schedule builder construction time:\s*([\d.]+) ms", eng)
            eval_ms  = grabf(r"GPU kernel:.*=\s*([\d.]+) ms/vector", eng)
            gram_s   = grabf(r">>\s*total time:\s*([\d.]+)", gr.get(key, "")) \
                       or grabf(r"Compression time.*?([\d.]+)\s*sec", gr.get(key, ""))
            if build_ms is None:
                print(f"  [warn] no +pt build time for {key}", file=sys.stderr)
                continue
            ratio = (build_ms / eval_ms) if eval_ms else None
            rows.append(
                f"{tex} & {f2(gram_s)} & {f2(build_ms)} & "
                f"{('%.0f' % ratio) if ratio else '--'}\\\\")
    write_table("tab_build.tex",
                [r"matrix & grammar RePair (s) & $+\text{pt}$ build (ms) & $\approx$ mat-vecs\\"],
                rows, "lrrr")


def gen_fig_space_dat():
    """fig:space data: nnz, engine analytic MB, cuSPARSE analytic MB, oom flag (0/1).

    Both series are the *analytic* device footprint so the two are compared on the same
    metric (as tab:geno does). Plotting the engine's measured peak against cuSPARSE's
    analytic figure would flatter the engine -- the peak runs well under the analytic
    bound on some datasets (e.g. geno20: 12.0 vs 25.2 MB).
    """
    se = _space_energy_map()
    if not se:
        return
    rows = []
    for key, tex, r, c, nnz in GENO + [CROSSOVER]:
        d = se.get(key)
        if not d:
            continue
        eng = (d["eng_mb"] or d["eng_pkmb"])   # analytic; peak only as a fallback
        cus = d["cus_mb"]
        oom = 1 if d["oom"] else 0
        # For the OOM point plot the analytic CSR estimate if present, else leave blank.
        cus_s = f"{cus:.1f}" if cus else "nan"
        rows.append(f"{(d['nnz'] or nnz):.0f} {eng:.1f} {cus_s} {oom} {key}")
    write_dat("fig_space.dat", "nnz eng_mb cus_mb oom label", rows)


# --------------------------------------------------------------------------------------
# Table 4 -- tab:geno_spmm  (best batched ms/vector, engine vs best cuSPARSE alg)
# Figure 4 -- fig:batched   (geno22full time vs B)
#   logs: manuscript/logs/geno_spmm.log         (engine sweep: "GPU batched (B=NN): ..., X ms/vector")
#         manuscript/logs/geno_cusparse_alg.log (alg sweep: "cuSPARSE SpMM (B=NN, alg): ..., X ms/vector")
# --------------------------------------------------------------------------------------
def _engine_spmm(blocks, key):
    """List of (B, ms/vector) from an engine sweep block."""
    out = []
    for m in re.finditer(r"GPU batched \(B=(\d+)\):.*?,\s*([\d.]+) ms/vector", blocks.get(key, "")):
        out.append((int(m.group(1)), float(m.group(2))))
    return out


ALG_ABBR = {"default": "def", "csr_alg2": "a2", "csr_alg3": "a3"}

def _cusparse_alg(blocks, key):
    """List of (B, alg, ms/vector) from a cuSPARSE alg-sweep block."""
    out = []
    for m in re.finditer(r"cuSPARSE SpMM \(B=(\d+),\s*(\w+)\):.*?,\s*([\d.]+) ms/vector", blocks.get(key, "")):
        out.append((int(m.group(1)), m.group(2), float(m.group(3))))
    return out


def gen_tab_geno_spmm():
    eng_b = read_blocks("geno_spmm.log")
    cus_b = read_blocks("geno_cusparse_alg.log")
    if not eng_b or not cus_b:
        return
    rows = []
    for group, hdr in ((GENO_REAL, REAL_HDR % 4), (GENO_SYNTH, SYNTH_HDR % 4)):
        rows.append(hdr)
        for key, tex, r, c, nnz in group:
            es = _engine_spmm(eng_b, key)
            cs = _cusparse_alg(cus_b, key)
            if not es or not cs:
                print(f"  [warn] no spmm for {key}", file=sys.stderr)
                continue
            eB, ems = min(es, key=lambda t: t[1])
            cB, calg, cms = min(cs, key=lambda t: t[2])
            ratio = ems / cms if cms else None
            e_str = (f"{ems:.3f}" if ems < 1 else f"{ems:.2f}")
            c_str = (f"{cms:.3f}" if cms < 1 else f"{cms:.2f}")
            rows.append(
                f"{tex} & {e_str} ({eB}) & {c_str} ({ALG_ABBR.get(calg, calg)},{cB}) & "
                f"{('%.1f' % ratio) if ratio else '--'}\\\\")
    write_table("tab_geno_spmm.tex",
                [r"matrix & engine ($B$) & cuSPARSE (alg,$B$) & $\times$cuS\\"],
                rows, "lrrr")


def gen_fig_batched_dat(key="geno22full"):
    """fig:batched data for one matrix: B, engine, cusparse best-of-default/alg2, cusparse alg3."""
    eng_b = read_blocks("geno_spmm.log")
    cus_b = read_blocks("geno_cusparse_alg.log")
    if not eng_b or not cus_b:
        return
    eng = dict(_engine_spmm(eng_b, key))
    cus = {}
    for B, alg, ms in _cusparse_alg(cus_b, key):
        cus.setdefault(B, {})[alg] = ms
    rows = []
    for B in sorted(eng):
        c = cus.get(B, {})
        # ALG_DEFAULT/ALG2 track together; report the better of the two.
        defalg2 = [c[a] for a in ("default", "csr_alg2") if a in c]
        alg3 = c.get("csr_alg3")
        da = min(defalg2) if defalg2 else float("nan")
        rows.append(f"{B} {eng[B]:.4f} {da:.4f} {alg3:.4f}" if alg3 else f"{B} {eng[B]:.4f} {da:.4f} nan")
    write_dat("fig_batched.dat", "B engine cusparse_defalg2 cusparse_alg3", rows)


# --------------------------------------------------------------------------------------
# Table 5 -- tab:graph_struct  (Wikidata structural figures)
#   log: manuscript/logs/graph_struct.log
# --------------------------------------------------------------------------------------
def gen_tab_graph_struct():
    blocks = read_blocks("graph_struct.log")
    if not blocks:
        return
    rows = []
    for key, tex, r, c, nnz in GRAPH:
        s = parse_struct(blocks.get(key, ""))
        if not s:
            print(f"  [warn] no struct for {key}", file=sys.stderr)
            continue
        # Single-column table: K/M-abbreviated numbers, no dimension suffix
        # (dimensions are given in the sec:semiring prose). nnz comes from the
        # GRAPH config -- the struct log does not emit an edge count.
        rows.append(f"{tex} & {abbr(nnz)} & {abbr(s['raw'])} & {abbr(s['total'])} & "
                    f"{s['L']} & {abbr(s['w'])} & {abbr(s['pt'])}\\\\")
    write_table("tab_graph_struct.tex",
                [r"relation & nnz & $|\mathcal{R}|$ & total & $L$ & $w^{*}$ & $+\text{pt}$\\"],
                rows, "lrrrrrr")


# --------------------------------------------------------------------------------------
# Table 6 -- tab:graph  (Boolean & Tropical right product across Wikidata relations)
#   log: manuscript/logs/graph_semiring.log
#        blocks keyed ## <key>_<semiring>  (e.g. wd_sports_team_boolean).
#        engine single ms:  "GPU kernel: ... = X ms/vector"
#        engine B=16 ms:    "GPU batched (B=16): ..., X ms/vector"
#        GraphBLAS ms:      "GraphBLAS ...: X ms/vector"   (or "GB mxv: X ms")
#        cuSPARSE (bool):   "Average time: X ms/vector"
#        analytic MB:       "MEM analytic bytes: N"  (engine); CSR analytic from cusparse
# --------------------------------------------------------------------------------------
def _graph_sem(blocks, key, sem):
    b = blocks.get(f"{key}_{sem}", "")
    if not b:
        return None
    return dict(
        eng_single = grabf(r"GPU kernel:.*=\s*([\d.]+) ms/vector", b),
        eng_b16    = grabf(r"GPU batched \(B=16\):.*?,\s*([\d.]+) ms/vector", b),
        # The manuscript's GraphBLAS baseline is the 20-thread run (see tab:graph
        # caption); the log now emits both nthreads=1 and nthreads=20, so prefer the
        # 20-thread line and fall back to the first GraphBLAS line for single-line logs.
        gb         = (grabf(r"(?:GraphBLAS|GB)[^\n]*?nthreads=\s*20[^\n]*?([\d.]+)\s*ms", b)
                      or grabf(r"(?:GraphBLAS|GB)[^\n]*?([\d.]+) ms(?:/vector)?", b)),
        cus        = grabf(r"Average time:\s*([\d.]+) ms/vector", b),
        eng_mb     = (grabi(r"MEM analytic bytes:\s*(\d+)", b) or 0)/1e6 or None,
    )


def gen_tab_graph():
    blocks = read_blocks("graph_semiring.log")
    struct = read_blocks("graph_struct.log") or {}
    if not blocks:
        return
    rows = []
    for key, tex, r, c, nnz in GRAPH:
        bo = _graph_sem(blocks, key, "boolean")
        tr = _graph_sem(blocks, key, "tropical")
        if not bo or not tr:
            print(f"  [warn] no semiring data for {key}", file=sys.stderr)
            continue
        s = parse_struct(struct.get(key, "")) or {}
        # cuSPARSE CSR analytic footprint: the boolean block has two "MEM analytic
        # bytes" lines (engine first, then cuSPARSE after the "Reconstructed CSR"
        # banner); take the cuSPARSE one.
        bo_raw = blocks.get(f"{key}_boolean", "")
        cus_half = re.split(r"Reconstructed CSR|--- CUSPARSE ---", bo_raw, maxsplit=1)
        cus_bytes = grabi(r"MEM analytic bytes:\s*(\d+)", cus_half[1]) if len(cus_half) > 1 else None
        cus_mb = (cus_bytes / 1e6) if cus_bytes else None
        rows.append(
            f"{tex} & {grp(nnz)} & {grp(s.get('raw'))} & "
            f"{f2(cus_mb)}/{f2(bo['eng_mb'])} & "
            f"{f3(bo['cus'])}/{f3(bo['eng_single'])}/{f3(bo['gb'])} & {f3(bo['eng_b16'])} & "
            f"{f3(tr['eng_single'])}/{f3(tr['gb'])} & {f3(tr['eng_b16'])} \\\\")
    header = [
        r"& & & analytic MB & \multicolumn{2}{c}{Boolean} & \multicolumn{2}{c}{Tropical}\\",
        r"\cmidrule(lr){5-6}\cmidrule(lr){7-8}",
        r"relation & nnz & $|R|$ & (CSR/eng) & single (CSR/eng/GB) & $B{=}16$ (eng) & single (eng/GB) & $B{=}16$ (eng) \\",
    ]
    write_table("tab_graph.tex", header, rows, "lrrr rr rr")


# --------------------------------------------------------------------------------------
# Table -- tab:graph_scale  (largest two Wikidata relations, bracketing compressibility)
#   log: manuscript/logs/graph_scale.log, one ## <key> block per relation, three sub-parts:
#     engine:   SEMIRING=boolean gpu_test -> struct (|R|,L), "MEM analytic bytes" (eng dev),
#               "GPU kernel: ... = X ms/vector" (Boolean engine time)
#     "--- CUSPARSE ---" cusparse_test -> "MEM analytic bytes" (CSR), "Average time" (Bool cuS)
#     "--- GRAMMAR ---"  matrepair -r -y --bool -> " >> REANS size: N bytes" (serialized)
# --------------------------------------------------------------------------------------
def gen_tab_graph_scale():
    blocks = read_blocks("graph_scale.log")
    if not blocks:
        return
    rows = []
    for key, tex, r, c, nnz in GRAPH_SCALE:
        b = blocks.get(key, "")
        eng, rest = (b.split("--- CUSPARSE ---", 1) + [""])[:2]
        cus, gram = (rest.split("--- GRAMMAR ---", 1) + [""])[:2]
        s = parse_struct(eng)
        if not s:
            print(f"  [warn] no struct for {key}", file=sys.stderr)
            continue
        eng_mb   = grabi(r"MEM analytic bytes:\s*(\d+)", eng)
        cus_mb   = grabi(r"MEM analytic bytes:\s*(\d+)", cus)
        bool_eng = grabf(r"GPU kernel:.*=\s*([\d.]+) ms/vector", eng)
        bool_cus = grabf(r"Average time:\s*([\d.]+) ms/vector", cus)
        reans_b  = grabi(r">>\s*REANS size:\s*(\d+)\s*bytes", gram)
        rows.append(
            f"{tex} & {grp(nnz)} & {grp(s['raw'])} & {s['L']} & "
            f"{mb(reans_b)} & {mb(cus_mb)} & {mb(eng_mb)} & "
            f"{f2(bool_eng)} / {f2(bool_cus)}\\\\")
    header = [
        r"relation & nnz & $|\mathcal{R}|$ & $L$ & REANS & CSR & eng.\ dev. & Bool eng/cuS\\",
        r" & & & & (MB) & (MB) & (MB) & (ms)\\",
    ]
    write_table("tab_graph_scale.tex", header, rows, "lrrrrrrr")


# --------------------------------------------------------------------------------------
def main():
    print("Extracting manuscript tables + figure data from manuscript/logs/ ...")
    generators = [
        gen_tab_geno_through, gen_tab_geno_time, gen_tab_geno, gen_tab_geno_spmm,
        gen_tab_build, gen_tab_graph_struct, gen_tab_graph, gen_tab_graph_scale,
        gen_fig_space_dat, gen_fig_batched_dat,
    ]
    for g in generators:
        try:
            g()
        except Exception as e:                      # keep going; one bad log != total failure
            print(f"  [ERR] {g.__name__}: {e}", file=sys.stderr)
    print("Done. (Missing logs were skipped; rerun after the corresponding experiment.)")


if __name__ == "__main__":
    main()
