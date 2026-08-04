#!/usr/bin/env python3
"""Generate the E4 bibliometric re-analysis notebook (plan P3-T1[b]).

Rewrite of the v1 bibliometric analysis. ALL percentages are computed from the
coded data (results/E4/bibliometric_coded.csv), never by hand-typed literals
(plan §1.5 root cause: cell 5 hand-summed nine 1.01% literals and dropped a
tenth). Correct denominators are used: N_screened = 100, N_analysed = 99,
N_non-arXiv = 88. Figures are regenerated with larger fonts (SoutoNeto §4) and
correct captions, and the publisher/OA stratification requested by the referee
is added with an explicit small-cell caveat.

Regenerate: python notebooks/_generators/build_e4_analysis_nb.py
"""
import json, os, textwrap

OUT = "/home/hugo_souto/Stuff/Submitted_Research/TISCA_paper/Test-Informed-Simulation-Count-Algorithm-TISCA/notebooks/E4_bibliometric_analysis.ipynb"
CELLS = []

def md(src):
    CELLS.append({"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)})

def code(src):
    CELLS.append({"cell_type": "code", "execution_count": None, "metadata": {},
                  "outputs": [], "source": src.splitlines(keepends=True)})

md(textwrap.dedent("""\
    # E4 - Bibliometric re-analysis (P3-T1[b])

    Recount and figures for the Section 2 bibliometric study, computed **from
    code** (results/E4/bibliometric_coded.csv). No percentage literal appears
    anywhere in this notebook: every number is emitted by a cell. Correct
    denominators: N_screened = 100, N_analysed = 99, N_non-arXiv = 88. The
    `justification_*` columns of the coded CSV are `REVIEW_REQUIRED` pending the
    per-paper reading pass (plan P3-T1[a]); figure code below is defensive about
    that.

    On Colab: drop this notebook in, run top to bottom. It offers `files.download`
    for the figures (standing rule). Runs on laptop or free CPU runtime.
    """))

code(textwrap.dedent("""\
    import os, csv, re
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt

    # locate the coded CSV (repo checkout preferred)
    cands = [
        "/content/Test-Informed-Simulation-Count-Algorithm-TISCA/results/E4/bibliometric_coded.csv",
        os.path.join(os.getcwd(), "results", "E4", "bibliometric_coded.csv"),
    ]
    path = next((p for p in cands if os.path.exists(p)), None)
    if path is None:
        raise SystemExit("bibliometric_coded.csv not found; run experiments/E4_bibliometrics/code_bibliometrics.py first")
    df = pd.read_csv(path, keep_default_na=False)

    # Larger fonts for print legibility (SoutoNeto 4)
    plt.rcParams.update({"font.size": 14, "axes.labelsize": 15,
                         "axes.titlesize": 16, "xtick.labelsize": 13,
                         "ytick.labelsize": 13, "figure.figsize": (8, 6)})

    N_S = len(df)                                   # screened
    dfA = df[df["J_report_reported"].astype(int) == 1].copy()  # analysed (numeric J)
    N_A = len(dfA)
    dfA["J"] = dfA["J_numeric"].astype(int)
    dfA["year"] = pd.to_numeric(dfA["year"], errors="coerce")
    dfA["is_arxiv_preprint"] = dfA["is_arxiv_preprint"].astype(int)
    dfNA = dfA[dfA["is_arxiv_preprint"] == 0]
    N_NA = len(dfNA)
    print(f"N_screened = {N_S}")
    print(f"N_analysed (numeric J) = {N_A}")
    print(f"N_non-arXiv = {N_NA}")
    print(f"arxiv preprints = {dfA['is_arxiv_preprint'].sum()}")
    """))

code(textwrap.dedent("""\
    # ---- Audited headline percentages, all from counts ----
    def pct(num, den): return 100.0 * num / den if den else float("nan")
    arxiv_s = int(dfA["is_arxiv_preprint"].sum())
    y = dfA["year"]
    print("arXiv share            :", round(pct(arxiv_s, N_A), 2), "%", f"({arxiv_s}/{N_A})")
    print("2021-2025              :", round(pct(int(((y >= 2021) & (y <= 2025)).sum()), N_A), 2), "%")
    print("2016-2025              :", round(pct(int(((y >= 2016) & (y <= 2025)).sum()), N_A), 2), "%")
    print("J == 1000              :", round(pct(int((dfA.J == 1000).sum()), N_A), 2), "%")
    print("J <= 500 (all)         :", round(pct(int((dfA.J <= 500).sum()), N_A), 2), "%")
    print("J <= 500 (non-arXiv)   :", round(pct(int((dfNA.J <= 500).sum()), N_NA), 2), "%")
    # every percentage above must equal REVISION_PLAN.md 1.5 (11.11/56.57/87.88/
    # 34.34/54.55/52.27)
    """))

md(textwrap.dedent("""\
    #### Figure 1 - J distribution (all analysed papers, N = 99)
    Computed from `value_counts` over the analysed subset; "Not mentioned" is
    excluded so the denominator is 99, not 100.
    """))

code(textwrap.dedent("""\
    vc = dfA["J"].value_counts(normalize=True) * 100
    vc = vc.sort_values(ascending=False)
    fig, ax = plt.subplots()
    ax.barh([str(int(x)) for x in vc.index], vc.values, color="#4C72B0")
    ax.invert_yaxis()
    ax.set_xlabel("Percentage of analysed studies (%)")
    ax.set_ylabel("Number of simulations J")
    ax.set_title("J distribution (N = %d)" % N_A)
    for i, v in enumerate(vc.values):
        ax.text(v + 0.4, i, f"{v:.1f}%", va="center", fontsize=12)
    fig.tight_layout()
    fig.savefig("Fig_bib_J.png", dpi=300)
    plt.show()
    """))

md(textwrap.dedent("""\
    #### Figure 2 - Publisher / venue distribution (N = 99)
    """))

code(textwrap.dedent("""\
    vp = dfA["publisher"].replace("", dfA["venue"]).replace("", "unlisted")
    vc2 = vp.value_counts(normalize=True) * 100
    vc2 = vc2.head(12).sort_values()
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(vc2.index, vc2.values, color="#C44E52")
    ax.set_xlabel("Percentage of analysed studies (%)")
    ax.set_ylabel("Publisher / venue")
    ax.set_title("Publisher distribution (N = %d)" % N_A)
    for i, v in enumerate(vc2.values):
        ax.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=11)
    fig.tight_layout(); fig.savefig("Fig_bib_publisher.png", dpi=300); plt.show()
    """))

md(textwrap.dedent("""\
    #### Figures 3a / 3b - Distribution of publication year
    """))

code(textwrap.dedent("""\
    years = dfA["year"].dropna().astype(int)
    fig, ax = plt.subplots(figsize=(10, 5))
    years.value_counts().sort_index().plot(kind="bar", ax=ax, color="#55A868")
    ax.set_xlabel("Publication year"); ax.set_ylabel("Number of analysed studies")
    ax.set_title("Year distribution (N = %d)" % len(years))
    fig.tight_layout(); fig.savefig("Fig_bib_year.png", dpi=300); plt.show()
    """))

code(textwrap.dedent("""\
    # Non-arXiv subset year distribution (the v1 caption wrongly said 89; it is 88)
    years_na = dfNA["year"].dropna().astype(int)
    fig, ax = plt.subplots(figsize=(10, 5))
    years_na.value_counts().sort_index().plot(kind="bar", ax=ax, color="#8172B3")
    ax.set_xlabel("Publication year"); ax.set_ylabel("Number (non-arXiv)")
    ax.set_title("Year distribution, non-arXiv (N = %d)" % len(years_na))
    fig.tight_layout(); fig.savefig("Fig_bib_year_nonarxiv.png", dpi=300); plt.show()
    """))

md(textwrap.dedent("""\
    #### Stratification by publisher type and the journal-OA split

    `@misc` (preprint) versus formally published, and - for the journal subset
    that the mechanical coder flags REVIEW_REQUIRED - the gold-OA / hybrid /
    subscription judgement is a journal-policy fact to be confirmed in the reading
    pass. The code below reports the mechanically derivable split (preprint vs
    conference vs journal-flagged) and leaves the OA subtype as an explicit
    REVIEW_REQUIRED count. Small-cell caveat: any stratum below ~5 papers is not
    interpreted.
    """))

code(textwrap.dedent("""\
    strata = dfA["publisher_type"].value_counts()
    print("Mechanically-derivable publisher-type split (analysed, N = %d):" % N_A)
    print(strata.to_string())
    journal_flagged = int((dfA["publisher_type"] == "journal-UNCLASSIFIED-OA").sum())
    print("\\nJournal entries awaiting the gold-OA/hybrid/subscription reading pass:",
          journal_flagged)
    print("Small-cell caveat: strata with < 5 papers are not interpreted.")
    """))

code(textwrap.dedent("""\
    # Colab download fallback per the standing rule:
    try:
        from google.colab import files
        for fn in ["Fig_bib_J.png", "Fig_bib_publisher.png", "Fig_bib_year.png",
                   "Fig_bib_year_nonarxiv.png"]:
            if os.path.exists(fn):
                files.download(fn)
        print("Downloaded E4 figures.")
    except Exception as e:
        print("(not on Colab / download skipped):", e)
    """))

nb = {"nbformat": 4, "nbformat_minor": 0,
      "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3"},
                   "language_info": {"name": "python"}},
      "cells": CELLS}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    json.dump(nb, f, indent=1)
for i, c in enumerate(CELLS):
    src = "".join(c["source"])
    assert not any(b in src for b in ["'''", '"""']), f"cell {i} triple-quote"
print("wrote", OUT, len(CELLS), "cells")
