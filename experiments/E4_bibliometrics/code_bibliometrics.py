#!/usr/bin/env python3
"""E4 bibliometric re-coding (plan P3-T1[a]).

Reads the original bibliometric spreadsheet (legacy/Bibliometric_Study/
Bibliometric_study.xlsx) and emits results/E4/bibliometric_coded.csv with the
P3-T1(a) schema. All columns that are DERIVABLE FROM THE SPREADSHEET ALONE are
populated programmatically; columns that require reading each paper (whether,
and why, the authors justify their choice of J) are deliberately left with an
explicit `REVIEW_REQUIRED` marker because that coding is a manual, per-paper
judgement call (the plan budgets ~7 h of reading, split across 5 reviewers of
20 papers each). This script therefore never fabricates justification data: it
produces the auditable mechanical base and a ready-to-fill manual shell.

Mechanical columns derived here:
  paper_id   : the BibTeX cite key
  bibtex     : verbatim BibTeX string
  entry_type : @article / @inproceedings / @misc / ...
  venue      : journal name if present, else booktitle
  publisher  : publisher field if present, else venue/publisher heuristic
  year       : publication year
  peer_reviewed : N for @misc (preprint, by construction); Y for any formally
                 published type (@article/@inproceedings/@incollection/@book)
  J_verbatim / J_numeric : the recorded "Number of simulations" (and conversion)
  publisher_type :
      "preprint"    for @misc
      "conference"  for @inproceedings / @incollection
      "journal-UNCLASSIFIED-OA"  for @article: the gold-OA/hybrid/subscription
                 split is a journal-policy fact the bibtex does not carry, so it
                 is flagged for manual assignment (SECTION requires reading).
  is_arxiv_preprint : 1 when entry_type == @misc (the paper's v1 arxiv share)

Rows whose `Number of simulations` is literally "Not mentioned" get
`J_report_reported = 0` and are treated as screened-but-not-analysed (plan §1.5
N_screened=100, N_analysed=99).

Manual coding fields (left REVIEW_REQUIRED, to be filled by the reading pass):
  J_coding_rule, n_scenarios, J_is_outer_replication, confounded_with,
  justification_given, justification_quote, justification_type,
  publisher_type_gold_OA / _hybrid / _subscription.

Usage:
  python experiments/E4_bibliometrics/code_bibliometrics.py
Output:
  results/E4/bibliometric_coded.csv
"""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(
    HERE, "..", "..", "legacy", "Bibliometric_Study", "Bibliometric_study.xlsx"
)
OUT = os.path.join(HERE, "..", "..", "results", "E4", "bibliometric_coded.csv")


def load_rows(path):
    import openpyxl
    wb = openpyxl.load_workbook(path)
    ws = wb["Sheet1"]
    data = list(ws.iter_rows(min_row=2, values_only=True))
    if not all(len(r) == 2 for r in data):
        raise ValueError("spreadsheet must have exactly two columns")
    return data


def parse_bibtex_fields(bibtex):
    # Both the entry-type and the cite-key patterns are anchored with `re.match`, so
    # a single leading newline in the spreadsheet cell silently emptied BOTH: rows 33
    # and 82 (`@InProceedings{pmlr-v202-wu23i}` and `{pmlr-v151-xu22c}`) came through
    # with no entry type and the fallback `row_NN` id, which then denied them the
    # `inproceedings` branch in `classify_venues.py` and left two PMLR proceedings
    # coded `hybrid/subscription`. The stored `bibtex` column was already stripped,
    # so the CSV looked correct and the defect was invisible downstream.
    bibtex = (bibtex or "").strip()

    def field(name):
        m = re.search(
            name + r"\s*=\s*(?:\{([^{}]*)\}|([\"`'][^\"`']*[\"`']|[A-Za-z0-9_.\-]+))",
            bibtex, re.IGNORECASE | re.DOTALL,
        )
        if m:
            return (m.group(1) or m.group(2)).strip()
        return None

    key = None
    m = re.match(r"@\w+\s*[({]\s*([^,\s]+)\s*,", bibtex, re.DOTALL)
    if m:
        key = m.group(1).strip()
    etype = re.match(r"@\s*(\w+)", bibtex)
    entry_type = etype.group(1).lower() if etype else ""
    year = field("year")
    journal = field("journal")
    booktitle = field("booktitle")
    publisher = field("publisher")
    venue = journal or booktitle or ""
    return {
        "key": key,
        "entry_type": entry_type,
        "year": year,
        "journal": journal or "",
        "booktitle": booktitle or "",
        "publisher": publisher or "",
        "venue": venue,
    }


VENUES = os.path.join(HERE, "..", "..", "results", "E4", "venue_openaccess.csv")


def load_venue_types(path=VENUES):
    """Per-venue publisher type from ``classify_venues.py`` (DOAJ-sourced).

    Returns ``{}`` when the table has not been built yet, in which case the
    journal rows stay ``journal-UNCLASSIFIED-OA`` exactly as before. The lookup is
    on the same normalised-title key the classifier writes.
    """
    if not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        return {r["venue_key"]: r for r in csv.DictReader(f)}


def venue_key(venue, publisher):
    def norm(t):
        t = (t or "").replace("\\&", "&").replace("&amp;", "&")
        t = re.sub(r"[^a-z0-9&]+", " ", t.lower())
        return re.sub(r"\s+", " ", t).strip()
    return norm(venue) or f"__novenue__{norm(publisher)}"


def classify(entry, venue_types=None):
    et = entry["entry_type"]
    hit = (venue_types or {}).get(venue_key(entry["venue"], entry["publisher"]))
    if hit and hit["publisher_type"] != "UNRESOLVED":
        # The DOAJ-backed table supersedes the entry-type heuristic wherever it has
        # a verdict: it also corrects proceedings published under an @article entry
        # type (AAAI) and separates book chapters from conference papers.
        #
        # `working-paper` is NOT asserted either way. An entry with no venue may be
        # an unrefereed series item (NBER t0283) or a refereed article whose BibTeX
        # simply lost its journal field (sloczynski2022doubly), and the mechanical
        # record cannot tell them apart. Inferring "Y" put an unrefereed NBER
        # technical working paper into the peer-reviewed count.
        peer = {"preprint": "N", "working-paper": "REVIEW_REQUIRED"}.get(
            hit["publisher_type"], "Y")
        return {"publisher_type": hit["publisher_type"],
                "publisher_family": hit["publisher_family"],
                "peer_reviewed": peer,
                "publisher_type_source": hit["source"],
                "doaj_listed": hit["doaj_listed"]}
    if et == "misc":
        return {"publisher_type": "preprint", "peer_reviewed": "N"}
    if et in ("inproceedings", "conference"):
        return {"publisher_type": "conference", "peer_reviewed": "Y"}
    if et in ("incollection", "book", "inbook"):
        return {"publisher_type": "book-chapter", "peer_reviewed": "Y"}
    if et == "article":
        return {"publisher_type": "journal-UNCLASSIFIED-OA", "peer_reviewed": "Y"}
    if et in ("techreport", "phdthesis"):
        return {"publisher_type": "working-paper", "peer_reviewed": "N"}
    return {"publisher_type": "unknown", "peer_reviewed": "REVIEW_REQUIRED"}


def to_numeric_j(j_verbatim):
    if j_verbatim is None:
        return None
    s = str(j_verbatim).strip()
    if not s or s.lower() == "not mentioned":
        return None
    try:
        return int(float(s))
    except ValueError:
        m = re.search(r"\d+", s.replace(",", ""))
        return int(m.group()) if m else None


def main():
    rows = load_rows(SRC)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    venue_types = load_venue_types()
    cols = [
        "paper_id", "bibtex", "entry_type", "venue", "publisher", "year",
        "peer_reviewed", "publisher_type", "publisher_family",
        "publisher_type_source", "doaj_listed", "is_open_access",
        "J_verbatim", "J_numeric", "J_report_reported", "J_coding_rule",
        "n_scenarios", "J_is_outer_replication", "confounded_with",
        "justification_given", "justification_quote", "justification_type",
        "is_arxiv_preprint",
    ]
    recs = []
    for bidx, (bibtex, j_verbatim) in enumerate(rows):
        entry = parse_bibtex_fields(bibtex or "")
        cls = classify(entry, venue_types)
        j_num = to_numeric_j(j_verbatim)
        arxiv = 1 if entry["entry_type"] == "misc" else 0
        recs.append({
            "paper_id": entry["key"] or f"row_{bidx + 1}",
            "bibtex": (bibtex or "").strip(),
            "entry_type": entry["entry_type"],
            "venue": entry["venue"],
            "publisher": entry["publisher"],
            "year": entry["year"],
            "peer_reviewed": cls["peer_reviewed"],
            "publisher_type": cls["publisher_type"],
            "publisher_family": cls.get("publisher_family", "REVIEW_REQUIRED"),
            "publisher_type_source": cls.get("publisher_type_source", "entry type only"),
            "doaj_listed": cls.get("doaj_listed", ""),
            "is_open_access": 1 if cls["publisher_type"] in ("gold-OA", "diamond-OA") else 0,
            "J_verbatim": "" if j_verbatim is None else str(j_verbatim),
            "J_numeric": "" if j_num is None else j_num,
            "J_report_reported": 1 if j_num is not None else 0,
            "J_coding_rule": "REVIEW_REQUIRED",
            "n_scenarios": "REVIEW_REQUIRED",
            "J_is_outer_replication": "REVIEW_REQUIRED",
            "confounded_with": "REVIEW_REQUIRED",
            "justification_given": "REVIEW_REQUIRED",
            "justification_quote": "",
            "justification_type": "REVIEW_REQUIRED",
            "is_arxiv_preprint": arxiv,
        })

    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(recs)

    n = len(recs)
    analysed = sum(1 for r in recs if r["J_report_reported"] == 1)
    not_reported = n - analysed
    jle500 = sum(1 for r in recs if r["J_numeric"] != "" and r["J_numeric"] <= 500)
    print(f"wrote {os.path.abspath(OUT)}")
    print(f"  screened          : {n}")
    print(f"  analysed (numeric J): {analysed}")
    print(f"  not mentioned / non-numeric: {not_reported}")
    print(f"  arxiv preprints   : {sum(1 for r in recs if r['is_arxiv_preprint'])}")
    print(f"  J <= 500 (n={analysed}): {jle500}  ({100*jle500/analysed:.2f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
