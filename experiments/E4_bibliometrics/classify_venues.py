#!/usr/bin/env python3
"""Assign a publisher type to every venue in the bibliometric sample (plan P3-T1[a]).

SoutoNeto's report asks for the distribution of J stratified by publisher type, and
by gold-OA versus hybrid/subscription in particular. That is a journal-policy fact
the BibTeX does not carry, so ``code_bibliometrics.py`` left it as
``journal-UNCLASSIFIED-OA``. Filling it in from memory would put ~40 unverifiable
claims about named journals into a published table, so the gold-OA determination is
instead taken from **DOAJ membership**, which is a citable, reader-checkable
criterion: DOAJ lists only fully open-access journals, and it records whether an APC
is charged (so gold and diamond can be separated).

Resulting `publisher_type`:

    preprint            @misc (arXiv)
    conference          @inproceedings/@incollection, plus proceedings published
                        under an @article entry type (AAAI), fixed here
    working-paper       NBER/JSTOR-style entries with neither journal nor
                        booktitle, whatever entry type the exporter chose
    gold-OA             journal listed in DOAJ, APC charged
    diamond-OA          journal listed in DOAJ, no APC
    hybrid/subscription journal not listed in DOAJ

The last category deliberately merges hybrid and pure subscription. Splitting them
is a per-journal, per-year policy question with no comparable public register, and
every commercial and society publisher in this sample runs some OA option, so a
split would be both unverifiable and uninformative. ``publisher_family`` carries the
orthogonal, fully mechanical split (commercial / society / university-press /
preprint-server / conference) that the stratification can also use.

DOAJ responses are cached in ``results/E4/doaj_cache.json`` so the classification is
reproducible without network access, and re-running never re-queries a resolved
venue. Titles must match exactly after normalisation; a fuzzy hit is recorded as
unmatched rather than guessed.

Usage:
  python experiments/E4_bibliometrics/classify_venues.py            # uses cache, queries misses
  python experiments/E4_bibliometrics/classify_venues.py --offline  # cache only, no network
Output:
  results/E4/venue_openaccess.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CODED = os.path.join(ROOT, "results", "E4", "bibliometric_coded.csv")
CACHE = os.path.join(ROOT, "results", "E4", "doaj_cache.json")
OUT = os.path.join(ROOT, "results", "E4", "venue_openaccess.csv")

DOAJ_API = "https://doaj.org/api/search/journals/"
USER_AGENT = "TISCA-revision/1.0 (bibliometric venue classification)"

# Publisher family is mechanical: it is read off the publisher string, not judged.
COMMERCIAL = ("Informa", "Taylor", "Wiley", "Elsevier", "Springer", "SAGE",
              "de Gruyter", "Frontiers", "Curran")
UNIVERSITY_PRESS = ("Oxford University Press", "Cambridge University Press",
                    "MIT Press", "University of Chicago Press")
SOCIETY = ("Institute of Mathematical Statistics", "American Statistical",
           "American Psychological Association", "American Educational Research",
           "Association for the Advancement of Artificial Intelligence",
           "Proceedings of the National Academy of Sciences", "Royal Statistical",
           "Institute of Statistical Science", "Bangladesh Academy of Sciences",
           "International Biometric Society")

# Venues whose entry_type in the source .bib misrepresents the venue kind. These are
# corrections of a mechanical parsing artefact, not editorial judgements.
ENTRY_TYPE_OVERRIDES = {
    "proceedings of the aaai conference on artificial intelligence": "conference",
    "advances in neural information processing systems": "conference",
}


def norm(title):
    """Lowercase, strip LaTeX escapes/punctuation, collapse whitespace."""
    t = (title or "").replace("\\&", "&").replace("&amp;", "&")
    t = re.sub(r"[^a-z0-9&]+", " ", t.lower())
    return re.sub(r"\s+", " ", t).strip()


def load_cache():
    if os.path.exists(CACHE):
        with open(CACHE, encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(cache):
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=1, sort_keys=True)


def query_doaj(title, timeout=25):
    """Return the DOAJ record whose title matches ``title`` exactly, else None."""
    url = DOAJ_API + urllib.parse.quote(title, safe="") + "?pageSize=20"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.load(resp)
    want = norm(title)
    for hit in payload.get("results", []):
        bib = hit.get("bibjson", {})
        candidates = [bib.get("title", "")] + list(bib.get("alternative_title", []) or [])
        if isinstance(bib.get("alternative_title"), str):
            candidates.append(bib["alternative_title"])
        if any(norm(c) == want for c in candidates if c):
            apc = bib.get("apc", {}) or {}
            return {
                "matched": True,
                "doaj_title": bib.get("title", ""),
                "issn": ";".join(filter(None, [bib.get("pissn"), bib.get("eissn")])),
                "apc_charged": bool(apc.get("has_apc")),
                "publisher": (bib.get("publisher", {}) or {}).get("name", ""),
                "oa_start": (bib.get("oa_start") or ""),
            }
    return {"matched": False, "n_results": payload.get("total", 0)}


def publisher_family(publisher, kind):
    if kind == "preprint":
        return "preprint-server"
    if kind == "conference":
        return "conference"
    p = publisher or ""
    for needle in UNIVERSITY_PRESS:
        if needle.lower() in p.lower():
            return "university-press"
    for needle in SOCIETY:
        if needle.lower() in p.lower():
            return "society"
    for needle in COMMERCIAL:
        if needle.lower() in p.lower():
            return "commercial"
    return "other"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--offline", action="store_true",
                    help="use the cache only; never query DOAJ")
    ap.add_argument("--coded", default=CODED)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args(argv)

    if not os.path.exists(args.coded):
        print(f"missing {args.coded}; run code_bibliometrics.py first", file=sys.stderr)
        return 1
    with open(args.coded, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    cache = load_cache()
    venues = {}
    for r in rows:
        et, venue, pub = r["entry_type"], (r["venue"] or "").strip(), (r["publisher"] or "").strip()
        key = norm(venue) or f"__novenue__{norm(pub)}"
        rec = venues.setdefault(key, {"venue": venue, "publisher": pub,
                                      "entry_types": set(), "n_papers": 0})
        rec["entry_types"].add(et)
        rec["n_papers"] += 1

    queried = 0
    out_rows = []
    for key, rec in sorted(venues.items()):
        venue, pub, ets = rec["venue"], rec["publisher"], rec["entry_types"]
        if "misc" in ets:
            kind, conf, src = "preprint", "high", "entry type @misc"
        elif not venue:
            # Tested BEFORE the book branch. With neither `journal` nor `booktitle`
            # there is nothing the entry can be a chapter *of*, so the entry type an
            # exporter happened to choose carries no information about the venue
            # kind. Abadie & Imbens (2002) is the case that forced this: it is NBER
            # technical working paper t0283, exported as `@book` with an
            # `institution` field and no venue, and the old ordering coded it
            # `book-chapter`. Genuine chapters (nie2021nonparametric, booktitle
            # "Handbook of Statistical Methods for Precision Medicine") keep a venue
            # and are unaffected.
            kind, conf, src = "working-paper", "high", "no journal or booktitle field"
        elif ets & {"incollection", "inbook", "book"}:
            kind, conf, src = "book-chapter", "high", "entry type @incollection/@book"
        elif ets & {"inproceedings"} or key in ENTRY_TYPE_OVERRIDES:
            kind, conf, src = "conference", "high", "entry type / proceedings title"
        else:
            hit = cache.get(key)
            if hit is None and not args.offline:
                try:
                    hit = query_doaj(venue)
                    queried += 1
                    time.sleep(1.0)          # be polite to a free public API
                except Exception as exc:     # noqa: BLE001 - recorded, not swallowed
                    hit = {"matched": None, "error": repr(exc)}
                cache[key] = hit
                save_cache(cache)
            if hit is None:
                kind, conf, src = "UNRESOLVED", "none", "not cached, offline mode"
            elif hit.get("matched") is None:
                kind, conf, src = "UNRESOLVED", "none", f"DOAJ query failed: {hit.get('error','')}"
            elif hit["matched"]:
                kind = "gold-OA" if hit.get("apc_charged") else "diamond-OA"
                conf, src = "high", f"DOAJ exact title match ({hit.get('issn','')})"
            else:
                kind, conf, src = "hybrid/subscription", "high", "not listed in DOAJ"

        h = cache.get(key) or {}
        out_rows.append({
            "venue_key": key,
            "venue": venue,
            "publisher": pub,
            "n_papers": rec["n_papers"],
            "entry_types": "|".join(sorted(ets)),
            "publisher_type": kind,
            "publisher_family": publisher_family(pub, kind),
            "doaj_listed": {"gold-OA": 1, "diamond-OA": 1}.get(kind, 0),
            "doaj_apc_charged": int(bool(h.get("apc_charged"))) if h.get("matched") else "",
            "doaj_issn": h.get("issn", "") if h.get("matched") else "",
            "confidence": conf,
            "source": src,
        })

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)

    counts = {}
    papers = {}
    for r in out_rows:
        counts[r["publisher_type"]] = counts.get(r["publisher_type"], 0) + 1
        papers[r["publisher_type"]] = papers.get(r["publisher_type"], 0) + r["n_papers"]
    print(f"wrote {os.path.abspath(args.out)}  ({len(out_rows)} venues, {queried} DOAJ queries)")
    for k in sorted(counts):
        print(f"  {k:22s} {counts[k]:3d} venues  {papers[k]:3d} papers")
    unresolved = [r for r in out_rows if r["publisher_type"] == "UNRESOLVED"]
    if unresolved:
        print(f"  UNRESOLVED venues needing a manual check: "
              f"{[r['venue'] for r in unresolved]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
