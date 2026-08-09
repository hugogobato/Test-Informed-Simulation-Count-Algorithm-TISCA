#!/usr/bin/env python3
"""Code and analyse replication-number justifications (revision plan P5.5-T2).

The mechanical bibliometric file is intentionally left unchanged.  This script
joins a row-level manual coding manifest to that file, writes the complete
row-level output, and derives all percentages and Wilson intervals from the
output.  The manual manifest is embedded below so that the coding decisions and
their source locations travel with the analysis script.  A stable corpus row,
rather than ``paper_id``, is used as the key because the original corpus has
duplicate cite keys.

The justification study is an availability-stratified census.  A paper is
eligible when its numeric-J record is either an arXiv preprint or a journal
article classified as gold or diamond open access in the existing DOAJ-backed
venue table.  Papers outside that stratum are retained in the output but are
not labelled unjustified.
"""

import csv
import math
import os


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
BASE = os.path.join(ROOT, "results", "E4", "bibliometric_coded.csv")
OUT = os.path.join(ROOT, "results", "E4", "bibliometric_justifications_coded.csv")
SUMMARY = os.path.join(ROOT, "results", "E4", "justification_summary.csv")


def coding(
    *,
    access_basis,
    replication_unit,
    count_scope,
    other_counts_reported,
    j_coding_rule,
    n_scenarios,
    outer_replication,
    confounded_with,
    justification_class,
    justification_criterion,
    source_cited,
    sensitivity_to_j,
    sensitivity_evidence,
    source_url,
    source_location,
    source_version_checked,
    quote,
    coding_notes,
    secondary_repetition="none reported",
):
    """Return one manually coded record using the P5.5-T2 vocabulary."""
    given = {
        "explicit": "yes",
        "implicit_or_convention": "yes",
        "unjustified": "no",
        "unclear_report": "unclear",
        "not_applicable": "not_applicable",
    }[justification_class]
    status = {
        "explicit": "justification_present",
        "implicit_or_convention": "justification_present",
        "unjustified": "missing_justification",
        "unclear_report": "unclear_reporting",
        "not_applicable": "not_applicable",
    }[justification_class]
    return {
        "access_basis": access_basis,
        "eligibility_status": "coded_accessible_numeric_J",
        "justification_eligible": "Y",
        "replication_unit": replication_unit,
        "stated_count": "",
        "reported_count_scope": count_scope,
        "other_counts_reported": other_counts_reported,
        "justification_class": justification_class,
        "justification_status": status,
        "justification_criterion": justification_criterion,
        "source_cited": source_cited,
        "sensitivity_to_J": sensitivity_to_j,
        "sensitivity_evidence": sensitivity_evidence,
        "source_url": source_url,
        "source_location": source_location,
        "source_version_checked": source_version_checked,
        "justification_given": given,
        "justification_quote": quote,
        "secondary_repetition": secondary_repetition,
        "coding_notes": coding_notes,
        "J_coding_rule": j_coding_rule,
        "n_scenarios": n_scenarios,
        "J_is_outer_replication": outer_replication,
        "confounded_with": confounded_with,
        "justification_type": justification_class,
        "double_coded": "N",
        "coding_reviewer": "HGS",
    }


ARXIV = "https://arxiv.org/abs/{}"


# Full-text coding for the 17 records in the pre-specified accessible stratum.
# The quotations identify the count report, not a claim that the count is
# statistically justified.  A design, DGP, prior, or tuning rationale is not
# counted as a replication-count justification unless it addresses J itself.
MANUAL = {
    4: coding(
        access_basis="publisher diamond-OA full text; arXiv version",
        replication_unit="outer simulated data sets, independent replications within each DGP",
        count_scope="per DGP/scenario",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="8 DGPs",
        outer_replication="Y",
        confounded_with="none reported",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 200-replication count was reported.",
        source_url="https://doi.org/10.1214/19-BA1195 | " + ARXIV.format("1706.09523"),
        source_location="Simulation protocol and Section 4.1",
        source_version_checked="published article and arXiv",
        quote="Results are based on 200 independent replications for each DGP.",
        coding_notes="The paper motivates the simulation design and DGPs, but gives no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale.",
    ),
    19: coding(
        access_basis="arXiv full text",
        replication_unit="outer simulated data sets",
        count_scope="per simulation setting",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="multiple appendix settings",
        outer_replication="Y",
        confounded_with="none reported",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-repetition count was reported.",
        source_url=ARXIV.format("2212.03578"),
        source_location="Appendix D, simulation description and tables",
        source_version_checked="arXiv version checked",
        quote="For all the analyses, we simulate 1,000 times a dataset of size n = 1,000.",
        coding_notes="The stated count is a repeated simulated-data count. No count-specific rationale is given.",
    ),
    20: coding(
        access_basis="arXiv full text",
        replication_unit="outer simulated data sets",
        count_scope="high-dimensional simulation setup; corpus J maps to this setup",
        other_counts_reported="500 simulations in the piecewise-polynomial setup",
        j_coding_rule="per_scenario",
        n_scenarios="multiple simulation setups",
        outer_replication="Y",
        confounded_with="independent test samples are a separate repeated quantity",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of either the 100- or 500-simulation count was reported.",
        source_url=ARXIV.format("2004.14497"),
        source_location="Section 4.4, Simulation Experiments",
        source_version_checked="arXiv version checked",
        quote="Figure 4b shows mean squared-errors ... across 100 simulations.",
        coding_notes="The paper reports more than one simulation count. The corpus value 100 is identifiable as the high-dimensional setup, while the paper also reports 500 in the baseline setup. Neither count receives a count-specific rationale.",
        secondary_repetition="500 independent test samples in the baseline setup",
    ),
    22: coding(
        access_basis="publisher gold-OA full text",
        replication_unit="outer simulated data sets",
        count_scope="per sample-size/scenario setting",
        other_counts_reported="5,000 posterior draws after burn-in within each replicate",
        j_coding_rule="per_scenario",
        n_scenarios="3 sample sizes",
        outer_replication="Y",
        confounded_with="posterior draws are an inner Monte Carlo quantity, not J",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100-replicate count was reported.",
        source_url="https://doi.org/10.3389/fams.2023.1122114",
        source_location="Simulation study, Methods and results",
        source_version_checked="publisher full text checked",
        quote="All results were summarized over N = 100 replicates.",
        coding_notes="The paper clearly separates 100 outer replicates from 5,000 posterior draws, but gives no rationale for the outer count.",
        secondary_repetition="5,000 posterior draws after burn-in",
    ),
    32: coding(
        access_basis="arXiv full text",
        replication_unit="outer simulated or semi-simulated data sets/runs",
        count_scope="per experiment",
        other_counts_reported="20 additional runs for parameter tuning",
        j_coding_rule="per_scenario",
        n_scenarios="multiple synthetic and semi-synthetic experiments",
        outer_replication="Y",
        confounded_with="random data splits and tuning runs are reported separately",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100 performance-repeat count was reported.",
        source_url=ARXIV.format("2202.01336"),
        source_location="Main experiments and Appendix",
        source_version_checked="arXiv version checked",
        quote="We report results based on 100 repeats.",
        coding_notes="The 20 additional runs are used for tuning and are not counted as sensitivity to J. No rationale for 100 is reported.",
        secondary_repetition="20 additional tuning runs",
    ),
    37: coding(
        access_basis="arXiv full text",
        replication_unit="simulation random seeds generating data and fitting runs",
        count_scope="per simulated dataset and experiment",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="3 simulated datasets plus real-data analysis",
        outer_replication="Y",
        confounded_with="seed also governs fitting randomness; the two roles are not separated",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the five-seed count was reported.",
        source_url=ARXIV.format("2407.05287"),
        source_location="Simulation results tables and Appendix E",
        source_version_checked="arXiv version checked",
        quote="Reported: Average RMSE ± standard deviation over 5 random seeds.",
        coding_notes="The paper reports random seeds rather than separately naming independent outer data sets and optimization seeds. The count is therefore retained as the reported outer simulation unit, with the seed role flagged as unresolved. No count rationale is given.",
    ),
    43: coding(
        access_basis="arXiv full text",
        replication_unit="outer simulated data sets",
        count_scope="per setting",
        other_counts_reported="independent test sample of size 5,000 per setting",
        j_coding_rule="per_scenario",
        n_scenarios="2 main settings plus additional simulations",
        outer_replication="Y",
        confounded_with="independent evaluation sample is a separate quantity",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100-repetition count was reported.",
        source_url=ARXIV.format("2410.06377"),
        source_location="Section 6, simulation study",
        source_version_checked="arXiv version checked",
        quote="The sample size for each setting was 500 and the simulation was repeated 100 times.",
        coding_notes="The independent test sample is not an outer replication count. No rationale for 100 is reported.",
        secondary_repetition="independent test sample of size 5,000",
    ),
    46: coding(
        access_basis="publisher gold-OA full text; arXiv version",
        replication_unit="outer simulated data sets",
        count_scope="per simulation setting",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="multiple settings and estimator comparisons",
        outer_replication="Y",
        confounded_with="test samples are separate from the simulated training data",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-simulation count was reported.",
        source_url="https://doi.org/10.1515/jci-2022-0070 | " + ARXIV.format("2210.08272"),
        source_location="Section 4.2 and Tables 1 and 2",
        source_version_checked="publisher full text and arXiv",
        quote="We evaluate the performance of each estimator across 1000 simulations on test samples of size 1000.",
        coding_notes="The paper reports a separate test sample size and does not give a count-specific rationale for 1,000.",
        secondary_repetition="test samples of size 1,000",
    ),
    49: coding(
        access_basis="arXiv full text",
        replication_unit="outer feature/outcome simulations conditional on a generated network",
        count_scope="main simulation setting; corpus J maps to the 80-run setting",
        other_counts_reported="40 repeated runs in a second setting",
        j_coding_rule="per_scenario",
        n_scenarios="2 main network-size settings",
        outer_replication="Y",
        confounded_with="same network G reused across runs; network-generation uncertainty is not resampled",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 80-run count was reported.",
        source_url=ARXIV.format("2410.11797"),
        source_location="Sections 5.2 and 5.3; Appendix C",
        source_version_checked="arXiv version checked",
        quote="The simulation consists of 80 repeated runs, all based on the same network G.",
        coding_notes="The repeated runs are not independent network draws because G is generated once and reused. The paper also reports 40 runs in another setting, but gives no count-specific rationale.",
    ),
    54: coding(
        access_basis="arXiv full text",
        replication_unit="outer simulated data sets",
        count_scope="per DGP",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="3 DGPs",
        outer_replication="Y",
        confounded_with="none reported",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100-run count was reported.",
        source_url=ARXIV.format("2203.10975"),
        source_location="Section 4 and Tables 1 and 3",
        source_version_checked="arXiv version checked",
        quote="PEHE and RMSE are averaged over 100 simulation runs.",
        coding_notes="No count-specific rationale is reported.",
    ),
    57: coding(
        access_basis="arXiv full text",
        replication_unit="outer replicate data sets",
        count_scope="main DGP; additional replicate sets used for design sensitivity",
        other_counts_reported="additional replicate data sets for prior and ICC sensitivity analyses",
        j_coding_rule="per_scenario",
        n_scenarios="main DGP plus prior/ICC sensitivity scenarios",
        outer_replication="Y",
        confounded_with="within-replicate model comparison uses the same data set for BCF and aBCF",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper varies prior and ICC settings, but does not vary or assess the 200-replication count.",
        source_url=ARXIV.format("2407.07067"),
        source_location="Section 4.1, DGP and simulation description",
        source_version_checked="arXiv version checked",
        quote="For our simulation study, we drew 200 replicate data sets using our DGP.",
        coding_notes="Prior and ICC sensitivity are design sensitivities, not sensitivity to J. No count-specific rationale is reported.",
    ),
    60: coding(
        access_basis="arXiv full text",
        replication_unit="outer DGP data sets",
        count_scope="per DGP/experiment",
        other_counts_reported="500 burn-in and 500 post-burn-in MCMC iterations within each Bayesian fit",
        j_coding_rule="per_scenario",
        n_scenarios="at least 2 simulation studies",
        outer_replication="Y",
        confounded_with="MCMC iterations are inner posterior-sampling quantities, not J",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-replication count was reported.",
        source_url=ARXIV.format("2407.11927"),
        source_location="Sections 3.2, 4.1, and 4.2",
        source_version_checked="arXiv version checked",
        quote="As before, we will run 1000 replications of the simulation study.",
        coding_notes="Posterior burn-in and sampling iterations are explicitly separated from the outer DGP replications. Visual convergence checks do not justify the outer count. No count-specific rationale is reported.",
        secondary_repetition="500 burn-in and 500 post-burn-in MCMC iterations",
    ),
    62: coding(
        access_basis="arXiv full text",
        replication_unit="outer collections of simulated data sets",
        count_scope="four main scenarios in the corpus version",
        other_counts_reported="current arXiv v3 reports 100 collections per scenario; corpus records J=50",
        j_coding_rule="unclear_version_mapping",
        n_scenarios="4 main scenarios, plus supplementary scenarios",
        outer_replication="Y",
        confounded_with="none reported",
        justification_class="unclear_report",
        justification_criterion="unclear",
        source_cited="unclear",
        sensitivity_to_j="unclear",
        sensitivity_evidence="The accessible version cannot establish whether the corpus value 50 was a prior-version count or a coding error.",
        source_url=ARXIV.format("2504.03480"),
        source_location="Section 4, Simulation Study, and Supplementary Appendix C",
        source_version_checked="arXiv v3, dated 2026-02-08",
        quote="For each of the four scenarios, we generated 100 collections of datasets.",
        coding_notes="The pre-specified corpus row records J=50, whereas the accessible arXiv v3 reports 100. The count-to-source mapping therefore requires author verification. This row is not treated as unjustified.",
    ),
    66: coding(
        access_basis="arXiv full text",
        replication_unit="outer matched simulated data sets",
        count_scope="per simulation setting",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="multiple matching and treatment-effect settings",
        outer_replication="Y",
        confounded_with="matching and balance filtering determine the accepted-data-set construction",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-data-set count was reported.",
        source_url=ARXIV.format("2308.02005"),
        source_location="Simulation sections and Appendix",
        source_version_checked="arXiv version checked",
        quote="We generate 1000 matched datasets.",
        coding_notes="The paper justifies the matching/balance criteria, not the number of accepted matched data sets. No count-specific rationale is reported.",
    ),
    70: coding(
        access_basis="publisher gold-OA full text; arXiv version",
        replication_unit="outer Monte Carlo trials",
        count_scope="per simulation experiment; corpus J maps to the primary 100,000-trial results",
        other_counts_reported="20,000 trials in several MSE figures; 10 repeated batches for standard errors in some tables",
        j_coding_rule="per_scenario",
        n_scenarios="survey, ATE, policy, and semi-synthetic experiments",
        outer_replication="Y",
        confounded_with="10 batches used for standard errors are a secondary repetition, not part of J",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper uses different counts across experiments but does not assess sensitivity of a fixed experiment to J.",
        source_url="https://doi.org/10.1515/jci-2022-0019 | " + ARXIV.format("2106.07695"),
        source_location="Sections 5.1 and Appendix B",
        source_version_checked="publisher full text and arXiv",
        quote="Each entry is the average threshold chosen over 100,000 trials and standard errors are over 10 replications.",
        coding_notes="The reported 100,000 is a primary outer-trial count, but the paper also uses 20,000 in other figures. No power, precision, pilot-variance, precedent, or computational-budget rationale for any count is stated.",
        secondary_repetition="10 batches for standard errors in selected results",
    ),
    71: coding(
        access_basis="publisher gold-OA full text; arXiv version",
        replication_unit="outer simulated clustered RCT data sets",
        count_scope="per specification",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="multiple cluster and service-receipt specifications",
        outer_replication="Y",
        confounded_with="clusters and individuals are within-data-set sample-size quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-simulation count was reported.",
        source_url="https://doi.org/10.1515/jci-2022-0033 | " + ARXIV.format("2205.02143"),
        source_location="Section 6, Simulations, and Appendix D",
        source_version_checked="publisher full text and arXiv",
        quote="We ran 1,000 simulations for each specification.",
        coding_notes="The paper motivates parameter choices and evaluates confidence-interval coverage, but gives no count-specific rationale for 1,000.",
    ),
    77: coding(
        access_basis="publisher gold-OA full text",
        replication_unit="outer simulated populations/data sets",
        count_scope="per simulation study",
        other_counts_reported="none reported",
        j_coding_rule="per_scenario",
        n_scenarios="4 settings",
        outer_replication="Y",
        confounded_with="clusters and units are within-data-set population quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-repetition count was reported.",
        source_url="https://doi.org/10.1515/jci-2022-0079",
        source_location="Section 4, Simulation, especially the simulation repetition paragraph",
        source_version_checked="publisher full text checked",
        quote="The simulation study is repeated 1,000 times.",
        coding_notes="The paper reports finite-sample settings and coverage results but no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale.",
    ),
}


NEW_COLUMNS = [
    "corpus_row",
    "access_basis",
    "eligibility_status",
    "justification_eligible",
    "replication_unit",
    "stated_count",
    "reported_count_scope",
    "other_counts_reported",
    "justification_class",
    "justification_status",
    "justification_criterion",
    "source_cited",
    "sensitivity_to_J",
    "sensitivity_evidence",
    "source_url",
    "source_location",
    "source_version_checked",
    "secondary_repetition",
    "coding_notes",
    "double_coded",
    "coding_reviewer",
]


def wilson(count, denominator, z=1.96):
    if denominator == 0:
        return "", ""
    p = count / denominator
    denominator_adj = 1 + z * z / denominator
    centre = (p + z * z / (2 * denominator)) / denominator_adj
    half = z * math.sqrt(
        p * (1 - p) / denominator + z * z / (4 * denominator * denominator)
    ) / denominator_adj
    return round(100 * max(0.0, centre - half), 4), round(100 * min(1.0, centre + half), 4)


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def inaccessible_defaults(row):
    numeric = bool(row["J_numeric"])
    if numeric:
        status = "not_coded_inaccessible"
        basis = "no OA/arXiv full text in the pre-specified accessible stratum"
    else:
        status = "not_applicable_J_not_reported"
        basis = "no numeric J reported"
    return {
        "access_basis": basis,
        "eligibility_status": status,
        "justification_eligible": "N",
        "replication_unit": "not coded",
        "stated_count": row["J_numeric"],
        "reported_count_scope": "not coded",
        "other_counts_reported": "not coded",
        "justification_class": "not_applicable",
        "justification_status": status,
        "justification_criterion": "not_applicable",
        "source_cited": "not_applicable",
        "sensitivity_to_J": "not_applicable",
        "sensitivity_evidence": "not coded",
        "source_url": "",
        "source_location": "",
        "source_version_checked": "",
        "justification_given": "not_applicable",
        "justification_quote": "",
        "secondary_repetition": "not coded",
        "coding_notes": "This row is retained for denominator auditability but was not read because it was outside the accessible stratum.",
        "J_coding_rule": row["J_coding_rule"],
        "n_scenarios": row["n_scenarios"],
        "J_is_outer_replication": row["J_is_outer_replication"],
        "confounded_with": row["confounded_with"],
        "justification_type": "not_applicable",
        "double_coded": "N",
        "coding_reviewer": "HGS",
    }


def build_rows():
    base_rows = read_csv(BASE)
    if len(base_rows) != 100:
        raise ValueError(f"expected 100 screened rows, found {len(base_rows)}")
    accessible_rows = {
        i
        for i, row in enumerate(base_rows, 1)
        if row["J_numeric"]
        and (
            row["is_arxiv_preprint"] == "1"
            or row["publisher_type"] in {"gold-OA", "diamond-OA"}
        )
    }
    if accessible_rows != set(MANUAL):
        raise ValueError(
            "manual coding keys do not match accessible rows: "
            f"accessible={sorted(accessible_rows)}, manual={sorted(MANUAL)}"
        )

    rows = []
    for corpus_row, base in enumerate(base_rows, 1):
        manual = MANUAL.get(corpus_row, inaccessible_defaults(base))
        if corpus_row in MANUAL:
            manual = dict(manual)
            manual["stated_count"] = base["J_numeric"]
        combined = {"corpus_row": corpus_row, **base, **manual}
        rows.append(combined)
    return rows


def write_rows(rows):
    base_fields = list(read_csv(BASE)[0].keys())
    fields = NEW_COLUMNS + [field for field in base_fields if field not in NEW_COLUMNS]
    with open(OUT, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows):
    eligible = [row for row in rows if row["justification_eligible"] == "Y"]
    summary = []

    for category in [
        "explicit",
        "implicit_or_convention",
        "unjustified",
        "unclear_report",
    ]:
        count = sum(row["justification_class"] == category for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census",
                "variable": "justification_class",
                "category": category,
                "count": count,
                "denominator": len(eligible),
                "percentage": round(100 * count / len(eligible), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    for category in ["yes", "no", "unclear"]:
        count = sum(row["justification_given"] == category for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census",
                "variable": "justification_given",
                "category": category,
                "count": count,
                "denominator": len(eligible),
                "percentage": round(100 * count / len(eligible), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    for category in ["yes", "no", "unclear"]:
        count = sum(row["sensitivity_to_J"] == category for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census",
                "variable": "sensitivity_to_J",
                "category": category,
                "count": count,
                "denominator": len(eligible),
                "percentage": round(100 * count / len(eligible), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    for status in [
        "coded_accessible_numeric_J",
        "not_coded_inaccessible",
        "not_applicable_J_not_reported",
    ]:
        count = sum(row["eligibility_status"] == status for row in rows)
        low, high = wilson(count, len(rows))
        summary.append(
            {
                "analysis_set": "all screened corpus",
                "variable": "eligibility_status",
                "category": status,
                "count": count,
                "denominator": len(rows),
                "percentage": round(100 * count / len(rows), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    fields = [
        "analysis_set",
        "variable",
        "category",
        "count",
        "denominator",
        "percentage",
        "wilson_95_low",
        "wilson_95_high",
    ]
    with open(SUMMARY, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)
    return summary


def main():
    rows = build_rows()
    write_rows(rows)
    summary = write_summary(rows)
    eligible = [row for row in rows if row["justification_eligible"] == "Y"]
    print(f"wrote {os.path.abspath(OUT)}")
    print(f"wrote {os.path.abspath(SUMMARY)}")
    print(f"screened={len(rows)} accessible_numeric_J={len(eligible)}")
    for row in summary:
        if row["analysis_set"] == "accessible numeric-J census" and row["variable"] == "justification_class":
            print(
                f"{row['category']}: {row['count']}/{row['denominator']} "
                f"({row['percentage']:.2f}%, Wilson 95% CI "
                f"{row['wilson_95_low']:.2f}% to {row['wilson_95_high']:.2f}%)"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
