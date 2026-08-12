# Local-PDF justification coding (P5.5-T2)

Row-level coding of the replication-count justification for the 65 corpus
records that were read from PDFs the authors held, rather than from the
open-access stratum. The open-access stratum's 34 records are coded inline in
`../code_justifications.py`; these three manifests carry the rest, and
`load_local_pdf_coding()` merges them. Together the two strata cover 99 of the
100 numeric-`J` records in the corpus.

| File | Content |
|---|---|
| `local_pdf_folder_{1,2,3}.csv` | one row per corpus record, in the `justification_codebook.md` vocabulary |
| `local_pdf_folder_{1,2,3}.md` | the reading notes: methods, per-row rationale, and the quotations the codes rest on |

`corpus_row` is the one-based row of `results/E4/bibliometric_coded.csv` and is
the merge key, because the mechanical corpus contains duplicate cite keys. The
merge refuses to run if any corpus row appears in more than one manifest or in
both strata.

The PDFs themselves are **not** committed. They are third-party copyrighted
articles; the manifests record the filename, the source version read, and the
locating quotation so that a reader can obtain the same source and check the
code. Classification rules, the permitted values of every field, and the
declared primary reporting rule (`unclear_report` counts as not justified) are
in `results/E4/justification_codebook.md`.

Fields these manifests do not carry, principally the `J` coding rule and the
scenario count, keep their mechanical values from `bibliometric_coded.csv`: the
strata are merged onto the mechanical base, never used to overwrite it.

No independent second coder was available, so `double_coded` is `N` on every
row. This is stated as a limitation in the manuscript rather than presented as
an agreement statistic.
