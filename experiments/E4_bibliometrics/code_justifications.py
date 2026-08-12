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
eligible when its numeric-J record is either an arXiv preprint, a journal
article classified as gold or diamond open access in the existing DOAJ-backed
venue table, or one of the manually verified conference full-text sources in
``FULLTEXT_ACCESSIBLE_ROWS``.  Papers outside that stratum are retained in the
output but are not labelled unjustified.
"""

import csv
import math
import os


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
BASE = os.path.join(ROOT, "results", "E4", "bibliometric_coded.csv")
OUT = os.path.join(ROOT, "results", "E4", "bibliometric_justifications_coded.csv")
SUMMARY = os.path.join(ROOT, "results", "E4", "justification_summary.csv")
MANUAL_CODING = os.path.join(HERE, "manual_coding")


#: Declared primary reporting rule: a record counts as justified only when the
#: accessible source states or cites a reason for the count. ``unclear_report`` is
#: therefore counted with ``unjustified``, not set aside, because the claim being
#: reported is "the paper gives the reader a reason for J", and a report the reader
#: cannot resolve does not give one. The rule is deliberately one-directional: it
#: can only understate how often a justification exists, so the headline is a
#: conservative lower bound on justification prevalence rather than a figure that
#: benefits from ambiguity. The four-way ``justification_class`` is retained in the
#: row-level file and reported disaggregated alongside every collapsed percentage,
#: so nothing is silently merged and the split can be undone by any reader.
PRIMARY_RULE = {
    "explicit": "justified",
    "implicit_or_convention": "justified",
    "unjustified": "not_justified",
    "unclear_report": "not_justified",
    "not_applicable": "not_applicable",
}

#: The coarse three-way field, kept for continuity with the original codebook. It
#: still separates ``unclear`` from ``no``; only ``justification_primary`` collapses
#: them, and only for reporting.
GIVEN_RULE = {
    "explicit": "yes",
    "implicit_or_convention": "yes",
    "unjustified": "no",
    "unclear_report": "unclear",
    "not_applicable": "not_applicable",
}


def normalise_class(justification_class):
    """Fold ``unclear_report`` into ``unjustified`` at the row level.

    A record whose count-to-source mapping cannot be resolved gives its reader no
    rationale for ``J``, which is the same reporting failure as a record that
    states none, so the two are one category in the released coding rather than
    two categories collapsed at summary time. The distinction is not discarded:
    the returned flag becomes the ``mapping_unresolved`` column, and the reading
    that produced it stays in ``coding_notes``, so the 15 affected records remain
    identifiable and the merge remains reversible.
    """
    if justification_class == "unclear_report":
        return "unjustified", "Y"
    return justification_class, "N" if justification_class != "not_applicable" else "not_applicable"


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
    justification_class, mapping_unresolved = normalise_class(justification_class)
    if mapping_unresolved == "Y":
        # The class now asserts that no rationale reaches the reader, so the
        # dependent fields cannot keep saying the rationale is unresolved.
        justification_criterion = ("none" if justification_criterion == "unclear"
                                   else justification_criterion)
        source_cited = "none" if source_cited == "unclear" else source_cited
    given = GIVEN_RULE[justification_class]
    status = {
        "explicit": "justification_present",
        "implicit_or_convention": "justification_present",
        "unjustified": "missing_justification",
        "unclear_report": "unclear_reporting",
        "not_applicable": "not_applicable",
    }[justification_class]
    return {
        "justification_primary": PRIMARY_RULE[justification_class],
        "access_basis": access_basis,
        "eligibility_status": "coded_accessible_numeric_J",
        "justification_eligible": "Y",
        "replication_unit": replication_unit,
        "stated_count": "",
        "reported_count_scope": count_scope,
        "other_counts_reported": other_counts_reported,
        "justification_class": justification_class,
        "mapping_unresolved": mapping_unresolved,
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


# Full-text sources in the corpus that were verified before coding.  These
# include proceedings pages or direct conference PDFs, plus arXiv versions of
# conference papers where the publisher page was not itself open. The
# OpenReview record in row 8 was also checked against an indexed author-hosted
# mirror. Row 35 is a subscription journal article whose full text the coder
# obtained and read after the automated access pass, so it is verified here
# rather than by the open-access venue rule.
FULLTEXT_ACCESSIBLE_ROWS = {
    5, 6, 7, 8, 9, 11, 15, 21, 25, 31, 33, 35, 42, 69, 75, 82, 89, 94
}


# Full-text coding for the 35 records in the expanded accessible stratum.
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
    5: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="outer simulated data sets or benchmark runs",
        count_scope="per simulation setting",
        other_counts_reported="500 independently generated test observations within each run; source reports 10 runs",
        j_coding_rule="unclear_count_mapping",
        n_scenarios="3 synthetic DGPs and multiple sample sizes",
        outer_replication="Y",
        confounded_with="test observations are within-run sample quantities",
        justification_class="unclear_report",
        justification_criterion="unclear",
        source_cited="none",
        sensitivity_to_j="unclear",
        sensitivity_evidence="The corpus records J=100, but the accessible paper states that results are averaged across 10 runs. The source does not establish what corpus value 100 refers to.",
        source_url="https://proceedings.mlr.press/v130/curth21a.html",
        source_location="Section 6.1, Synthetic Experiments",
        source_version_checked="PMLR proceedings PDF",
        quote="In all simulations, we evaluate performance on 500 independently generated test-observations, and average across 10 runs.",
        coding_notes="The accessible count statement does not match the corpus J=100. This row is not treated as unjustified until the corpus-to-source mapping is resolved.",
    ),
    6: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="semi-synthetic IHDP outcome realizations",
        count_scope="IHDP benchmark",
        other_counts_reported="separate Jobs benchmark; train, validation, and test split proportions 63/27/10",
        j_coding_rule="per_benchmark",
        n_scenarios="IHDP and Jobs benchmarks",
        outer_replication="Y",
        confounded_with="data splits are within-realization quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No sensitivity analysis varying the 1,000 IHDP realizations was reported.",
        source_url="https://proceedings.mlr.press/v70/shalit17a.html",
        source_location="Section 5.1, Simulated outcome: IHDP",
        source_version_checked="PMLR proceedings PDF",
        quote="We average over 1000 realizations of the outcomes with 63/27/10 train/validation/test splits.",
        coding_notes="The paper reports the benchmark count but does not give a count-specific power, precision, pilot-variance, precedent, or computational-budget rationale.",
    ),
    7: coding(
        access_basis="NeurIPS proceedings direct PDF",
        replication_unit="semi-synthetic IHDP outcome realizations",
        count_scope="IHDP benchmark",
        other_counts_reported="101 retained ACIC 2018 data sets; 25 reruns for each ACIC estimation procedure",
        j_coding_rule="per_benchmark",
        n_scenarios="IHDP and ACIC 2018 benchmarks",
        outer_replication="Y",
        confounded_with="ACIC data sets and reruns are separate benchmark quantities",
        justification_class="implicit_or_convention",
        justification_criterion="precedent",
        source_cited="Shalit et al. (2017), cited as [SJS16]",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper contrasts benchmark settings but does not vary the 1,000-realization IHDP count as a sensitivity analysis.",
        source_url="https://proceedings.neurips.cc/paper_files/paper/2019/file/8fb5f8be2aa9d6c64a04e3ab9f63feee-Paper.pdf",
        source_location="Section 5.1, Setup",
        source_version_checked="NeurIPS 2019 proceedings PDF",
        quote="Following [SJS16], we use 1000 realizations from the NPCI package.",
        coding_notes="This is a cited benchmark precedent rather than a new precision or power calculation for J.",
        secondary_repetition="25 estimation reruns for each ACIC procedure",
    ),
    8: coding(
        access_basis="OpenReview full-text record; indexed author-hosted mirror",
        replication_unit="semi-synthetic benchmark runs, exact outer unit not fully established",
        count_scope="benchmark experiments",
        other_counts_reported="IHDP benchmark and multiple data-generating settings",
        j_coding_rule="unclear_count_mapping",
        n_scenarios="multiple benchmark settings",
        outer_replication="REVIEW_REQUIRED",
        confounded_with="not established from the accessible indexed copy",
        justification_class="unclear_report",
        justification_criterion="unclear",
        source_cited="unclear",
        sensitivity_to_j="unclear",
        sensitivity_evidence="The corpus records J=100, but the OpenReview endpoint was rate-limited during retrieval and the indexed copy did not expose a reliable count-to-source mapping for independent outer replications.",
        source_url="https://openreview.net/forum?id=HkxBJT4YvB | https://papersdb.cs.ualberta.ca/~papersdb/uploaded_files/1194/paper_CausalML_NeurIPS_2019.pdf",
        source_location="Empirical evaluation sections and benchmark description",
        source_version_checked="OpenReview metadata and indexed mirror",
        quote="No reliable count-specific statement was located in the accessible indexed copy.",
        coding_notes="The record is retained as readable-access evidence but is not collapsed into unjustified because the corpus J mapping remains unresolved.",
    ),
    9: coding(
        access_basis="NeurIPS Datasets and Benchmarks direct PDF",
        replication_unit="IHDP semi-synthetic outcome realizations",
        count_scope="IHDP-100 benchmark",
        other_counts_reported="IHDP-100 and IHDP-1000 are discussed; models are also repeated 5 times with different seeds within each run",
        j_coding_rule="per_benchmark",
        n_scenarios="IHDP case study and modified IHDP settings",
        outer_replication="Y",
        confounded_with="model-seed repetitions are secondary to the benchmark realizations",
        justification_class="implicit_or_convention",
        justification_criterion="precedent",
        source_cited="Curth et al. (2021) cite Shalit et al. (2017) as the IHDP benchmark source",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper discusses IHDP-100 and IHDP-1000 benchmark conventions but does not re-run its analysis under alternative J values.",
        source_url="https://datasets-benchmarks-proceedings.neurips.cc/paper_files/paper/2021/file/2a79ea27c279e471f4d180b08d62b00a-Paper-round2.pdf",
        source_location="Section 3, IHDP benchmarking practice, and Section 3 experimental setup",
        source_version_checked="NeurIPS Datasets and Benchmarks 2021 PDF",
        quote="We use [2]'s IHDP-100 dataset (100 realizations of the DGP) and report out-of-sample performance.",
        coding_notes="The selected count is inherited from a named benchmark. The paper also examines how results vary across the 100 realizations, but does not present a direct precision calculation for choosing 100.",
        secondary_repetition="Five model runs with different seeds per benchmark realization",
    ),
    11: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="outer simulated data sets",
        count_scope="Simulation 1 at n=500",
        other_counts_reported="50 replications for large-data simulations; MCMC iterations and burn-in are inner fitting quantities",
        j_coding_rule="per_scenario",
        n_scenarios="3 simulation studies",
        outer_replication="Y",
        confounded_with="MCMC iterations are inner posterior-sampling quantities, not J",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper reports 50 and 200 replications in different simulation studies but does not vary a fixed study's J to assess Monte Carlo sensitivity.",
        source_url="https://proceedings.mlr.press/v206/krantsevich23a.html",
        source_location="Section 5, Simulation Studies",
        source_version_checked="PMLR proceedings PDF",
        quote="For each of the methods, we averaged the results on the three metrics over 200 independent replications.",
        coding_notes="The paper reports different replication counts across studies and gives no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale.",
    ),
    15: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="benchmark evaluation runs",
        count_scope="per benchmark dataset",
        other_counts_reported="hyperparameter tuning, training iterations, and batch sizes are separate quantities",
        j_coding_rule="per_benchmark",
        n_scenarios="IHDP, Twins, and Jobs benchmarks",
        outer_replication="Y",
        confounded_with="training iterations and tuning repetitions are inner or calibration quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 10 benchmark runs was reported.",
        source_url="https://proceedings.mlr.press/v206/kuzmanovic23a.html",
        source_location="Section 5, Experiments, and Table 1 caption",
        source_version_checked="PMLR proceedings PDF",
        quote="Table 1: Results of experiments on three benchmark datasets (mean averaged over 10 runs ± standard deviation).",
        coding_notes="The paper reports the repeated-run count but gives no count-specific rationale.",
    ),
    21: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="not established; the corpus J=100 maps to an evaluation grid rather than outer data sets",
        count_scope="100 regular evaluation intervals in the synthetic metric calculation",
        other_counts_reported="NOuter=5000 for Monte Carlo inference; 20 objects with 10 instances in synthetic data; 25 samples per state in NEEC",
        j_coding_rule="not_outer_replication",
        n_scenarios="synthetic, IHDP, and NEEC benchmarks",
        outer_replication="N",
        confounded_with="evaluation grid, inference draws, and within-data-set object/instance counts",
        justification_class="unclear_report",
        justification_criterion="unclear",
        source_cited="none",
        sensitivity_to_j="unclear",
        sensitivity_evidence="The accessible paper reports 100 regular evaluation intervals, not 100 independent outer simulation replications. The corpus J label is therefore not established.",
        source_url="https://proceedings.mlr.press/v119/witty20a.html",
        source_location="Section 6, Experiments, and metric definition",
        source_version_checked="PMLR proceedings PDF",
        quote="We average over 100 regular intervals between the 5th and 95th percentile of treatment assignment in the observational data.",
        coding_notes="This is a corpus coding error or unresolved mapping, not evidence that the paper lacked a replication-count justification.",
    ),
    25: coding(
        access_basis="AAAI proceedings full text; arXiv version",
        replication_unit="outer synthetic or semi-synthetic training-data experiments",
        count_scope="per experimental setup",
        other_counts_reported="100 experiments; training sample sizes n=200 or n=500; test size m=100; ACIC source has 1000 observations",
        j_coding_rule="unclear_count_mapping",
        n_scenarios="four binary-treatment synthetic setups, four continuous-treatment setups, and ACIC semi-synthetic data",
        outer_replication="Y",
        confounded_with="corpus J=500 maps to training sample size n=500, not the outer experiment count",
        justification_class="unclear_report",
        justification_criterion="unclear",
        source_cited="Nie and Wager (2021) for synthetic setups; Shimoni et al. (2018) for ACIC",
        sensitivity_to_j="unclear",
        sensitivity_evidence="The accessible paper reports 100 experiments, but the corpus J=500 corresponds to one of the training sample sizes. The source therefore does not support a clean mapping from corpus J to the outer replication count.",
        source_url="https://ojs.aaai.org/index.php/AAAI/article/download/30025/31803 | " + ARXIV.format("2312.10435"),
        source_location="Section 6, Experiments, and Appendix F",
        source_version_checked="AAAI-24 proceedings PDF and arXiv version",
        quote="We conducted 100 experiments ... In all experiments, we used n = 200 or n = 500 observations as training data.",
        coding_notes="The corpus extraction appears to have captured the n=500 training-sample size rather than the 100 outer experiments. This is a mapping problem, not evidence that the reported replication count lacks a rationale.",
    ),
    31: coding(
        access_basis="NeurIPS proceedings direct PDF",
        replication_unit="IHDP semi-synthetic benchmark realizations and simulation runs",
        count_scope="IHDP-100 setup and selected 100-run result summaries",
        other_counts_reported="101 main simulation settings; many plots average across 10 runs",
        j_coding_rule="per_benchmark",
        n_scenarios="Setups A-D, IHDP, ACIC2016, and Twins",
        outer_replication="Y",
        confounded_with="10-run summaries and 101 simulation settings are separate design quantities",
        justification_class="implicit_or_convention",
        justification_criterion="precedent",
        source_cited="Shalit et al. (2017), cited as the IHDP-100 benchmark source",
        sensitivity_to_j="no",
        sensitivity_evidence="The paper reports 10-run and 100-run summaries but does not vary the IHDP-100 count as a direct sensitivity analysis.",
        source_url="https://proceedings.neurips.cc/paper_files/paper/2021/file/8526e0962a844e4a2f158d831d5fddf7-Paper.pdf",
        source_location="Section 5.1, Experimental setup, and Figures 3-8",
        source_version_checked="NeurIPS 2021 proceedings PDF",
        quote="Here, we use the 90/10 train-test splits of [4]'s IHDP-100 benchmark to evaluate in- and out-of-sample performance.",
        coding_notes="The count is inherited from a named benchmark convention. The paper does not give a new precision or power calculation for 100.",
    ),
    33: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="independent synthetic and semi-synthetic data replications",
        count_scope="per synthetic or semi-synthetic dataset",
        other_counts_reported="10,000 training epochs and batch size 100 are inner optimization quantities",
        j_coding_rule="per_scenario",
        n_scenarios="synthetic and semi-synthetic PM-CMR experiments",
        outer_replication="Y",
        confounded_with="training epochs and batch size are inner optimization quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 10 independent replications was reported.",
        source_url="https://proceedings.mlr.press/v202/wu23i.html",
        source_location="Section 4, Experiments, and Appendix B.3",
        source_version_checked="PMLR proceedings PDF",
        quote="For the PM-CMR datasets, we conduct 10 independent replications.",
        coding_notes="The replication count is stated but no count-specific rationale is provided.",
    ),
    35: coding(
        access_basis="publisher full text (JEBS/AERA) verified and read by the coder",
        replication_unit="outer simulated data sets",
        count_scope="per simulation condition",
        other_counts_reported="not recorded at the time of the manual check",
        j_coding_rule="per_scenario",
        n_scenarios="not recorded at the time of the manual check",
        outer_replication="Y",
        confounded_with="none recorded",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 200-replication count was reported.",
        source_url="https://doi.org/10.3102/10769986221115446",
        source_location="Simulation study section",
        source_version_checked="published article",
        quote="",
        coding_notes="Read manually after the automated access pass; the corpus J=200 was confirmed against the published simulation study and no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale is given. The count report was confirmed but no verbatim quotation was transcribed.",
    ),
    42: coding(
        access_basis="NeurIPS proceedings direct PDF",
        replication_unit="outer simulated data sets",
        count_scope="per simulation setting",
        other_counts_reported="sample sizes and covariate dimensions vary across settings",
        j_coding_rule="per_scenario",
        n_scenarios="multiple theoretical and simulation settings",
        outer_replication="Y",
        confounded_with="sample size and covariate dimension are within-data-set design quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100-simulation count was reported.",
        source_url="https://proceedings.neurips.cc/paper_files/paper/2020/file/f75b757d3459c3e93e98ddab7b903938-Paper.pdf",
        source_location="Simulation section and empirical comparison",
        source_version_checked="NeurIPS 2020 proceedings PDF",
        quote="Estimation is measured via the root mean squared error (RMSE) averaged over 100 simulations.",
        coding_notes="The paper reports the count and estimator comparisons but gives no count-specific rationale.",
    ),
    69: coding(
        access_basis="AAAI proceedings full text; author-hosted PDF",
        replication_unit="independent synthetic data experiments",
        count_scope="per synthetic-data setting",
        other_counts_reported="sample sizes n=50, 100, 200; m=1000 or 5000 observations for variable-decomposition evaluation",
        j_coding_rule="per_scenario",
        n_scenarios="synthetic treatment-effect settings and real online-advertising application",
        outer_replication="Y",
        confounded_with="n and m are within-data-set sample-size quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 50 independent experiments was reported; changing n and m is not a sensitivity analysis for J.",
        source_url="https://ojs.aaai.org/index.php/AAAI/article/download/10480/10339 | https://pengcui.thumedialab.com/papers/ATE_DVD.pdf",
        source_location="ATE Estimation and Variables Decomposition subsections; Tables 1 and 2",
        source_version_checked="AAAI proceedings PDF and author-hosted PDF",
        quote="To evaluate the performance of our proposed method, we carry out the experiments 50 times independently.",
        coding_notes="The paper reports the outer count and uses it to compute Bias, SD, MAE, RMSE, TPR, and TNR, but gives no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale.",
    ),
    75: coding(
        access_basis="NeurIPS proceedings direct PDF",
        replication_unit="outer simulated data sets",
        count_scope="per simulation setting",
        other_counts_reported="10 iterations for train/test splits in the real-data application; optimization iterations are inner quantities",
        j_coding_rule="per_scenario",
        n_scenarios="4 simulation settings plus additional appendix settings",
        outer_replication="Y",
        confounded_with="test-set size and optimization iterations are separate quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 100 independent simulation replications was reported.",
        source_url="https://proceedings.neurips.cc/paper_files/paper/2024/file/0fd3d8093b3ba73d19b393a1326fdba7-Paper-Conference.pdf",
        source_location="Section 4, Simulation Studies, and Appendix I",
        source_version_checked="NeurIPS 2024 proceedings PDF",
        quote="Each simulation setting is replicated independently 100 times.",
        coding_notes="The paper provides extensive simulation details but no count-specific power, precision, pilot-variance, precedent, or computational-budget rationale for 100.",
    ),
    82: coding(
        access_basis="PMLR proceedings page and downloadable PDF",
        replication_unit="outer simulation replicates",
        count_scope="per simulation scenario",
        other_counts_reported="1000 bootstrap resamples for confidence intervals; sample sizes and bin counts vary",
        j_coding_rule="per_scenario",
        n_scenarios="RCT and observational settings with high-dimensional variants",
        outer_replication="Y",
        confounded_with="bootstrap resamples are inner uncertainty-estimation quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000 simulation replicates was reported.",
        source_url="https://proceedings.mlr.press/v151/xu22c.html",
        source_location="Section 3.1, Simulation Schema, and Section 4.1",
        source_version_checked="PMLR proceedings PDF",
        quote="We generate 1000 replicates for each simulation scenario.",
        coding_notes="The paper distinguishes simulation replicates from bootstrap resampling but gives no count-specific rationale for 1,000.",
    ),
    89: coding(
        access_basis="NeurIPS proceedings direct PDF",
        replication_unit="outer simulated outcome data sets",
        count_scope="IHDP simulated-outcome benchmark",
        other_counts_reported="500 repetitions for the synthetic data-generation experiment; parameter and dimension grids are design quantities",
        j_coding_rule="per_benchmark",
        n_scenarios="synthetic and IHDP benchmark experiments",
        outer_replication="Y",
        confounded_with="500 synthetic-data repetitions are a separate simulation study",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 200 IHDP outcome replications was reported.",
        source_url="https://proceedings.neurips.cc/paper_files/paper/2017/file/b2eeb7362ef83deff5c7813a67e14f0a-Paper.pdf",
        source_location="IHDP Dataset with Simulated Outcomes",
        source_version_checked="NeurIPS 2017 proceedings PDF",
        quote="We repeat such procedures for 200 times and generate 200 sets of simulated outcomes.",
        coding_notes="The paper reports the benchmark replication count and a separate 500-repetition synthetic experiment, but gives no count-specific rationale for 200.",
    ),
    94: coding(
        access_basis="arXiv full text; CIKM conference paper metadata",
        replication_unit="independent IHDP or synthetic data sets",
        count_scope="per dataset and simulation setting",
        other_counts_reported="synthetic data have 1000 units; 1000 observations are split into 750 control and 250 treatment units within each generated data set",
        j_coding_rule="per_scenario",
        n_scenarios="IHDP, modified IHDP, and synthetic simulation datasets",
        outer_replication="Y",
        confounded_with="unit-level sample sizes and treatment/control composition are within-data-set quantities",
        justification_class="unjustified",
        justification_criterion="none",
        source_cited="none",
        sensitivity_to_j="no",
        sensitivity_evidence="No variation of the 1,000-data-set count was reported; the sensitivity analysis varies model hyperparameters instead.",
        source_url="https://arxiv.org/abs/2009.06828 | https://doi.org/10.1145/3340531.3412037",
        source_location="Section 4.1, Datasets; Section 4.4, Results; Section 4.6, Sensitivity Analysis",
        source_version_checked="arXiv v2 PDF and CIKM 2020 metadata",
        quote="We repeat these procedures 1000 times to conduct evaluations of the uncertainty of estimates.",
        coding_notes="The paper reports 1,000 independent repetitions and separately states that the synthetic data set contains 1,000 units. The phrase about robust estimation describes the purpose of repeating, but does not supply a count-specific precision, power, pilot-variance, precedent, or computational-budget argument.",
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
    "mapping_unresolved",
    "justification_primary",
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
        basis = "no verified full text in the expanded accessible stratum"
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
        "mapping_unresolved": "not_applicable",
        "justification_primary": PRIMARY_RULE["not_applicable"],
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
        "coding_notes": "This row is retained for denominator auditability but was not read because it was outside the expanded accessible stratum.",
        "J_coding_rule": row["J_coding_rule"],
        "n_scenarios": row["n_scenarios"],
        "J_is_outer_replication": row["J_is_outer_replication"],
        "confounded_with": row["confounded_with"],
        "justification_type": "not_applicable",
        "double_coded": "N",
        "coding_reviewer": "HGS",
    }


def load_local_pdf_coding(base_rows):
    """Read the local-PDF stratum coding from ``manual_coding/local_pdf_folder_*.csv``.

    The open-access stratum is embedded in ``MANUAL`` above because it was coded
    first and its 34 records fit comfortably in the script. The remaining 65 records
    were read from PDFs the authors already held, and their row-level coding is kept
    as three CSV manifests rather than as another 800 lines of literals. The
    manifests use the same controlled vocabulary as ``coding()``; anything they do
    not recode (the ``J`` coding rule, the scenario count) stays at the mechanical
    value from ``bibliometric_coded.csv``, so the two strata are merged, never
    overwritten. Each folder's companion ``.md`` carries the reading notes.
    """
    coded = {}
    for folder in (1, 2, 3):
        path = os.path.join(MANUAL_CODING, f"local_pdf_folder_{folder}.csv")
        for record in read_csv(path):
            corpus_row = int(record["corpus_row"])
            if corpus_row in coded:
                raise ValueError(f"corpus row {corpus_row} coded twice in manifests")
            if not 1 <= corpus_row <= len(base_rows):
                raise ValueError(f"corpus row {corpus_row} is outside the corpus")
            filename = record.get("local_pdf_filename") or record["local_PDF_filename"]
            justification_class = record["justification_class"]
            if justification_class not in PRIMARY_RULE:
                raise ValueError(
                    f"row {corpus_row}: unknown class {justification_class!r}")
            justification_class, mapping_unresolved = normalise_class(
                justification_class)
            applicable = justification_class != "not_applicable"
            coded[corpus_row] = {
                "access_basis": f"local PDF: Remaining_Papers/Folder_{folder}/{filename}",
                "eligibility_status": ("coded_accessible_numeric_J" if applicable
                                       else "not_applicable_J_not_reported"),
                "justification_eligible": "Y" if applicable else "N",
                "replication_unit": record["replication_unit"],
                "stated_count": record["J_numeric"],
                "reported_count_scope": record["reported_count_scope"],
                "other_counts_reported": record["other_counts_reported"],
                "justification_class": justification_class,
                "mapping_unresolved": mapping_unresolved,
                "justification_primary": PRIMARY_RULE[justification_class],
                "justification_status": ("missing_justification"
                                         if mapping_unresolved == "Y"
                                         else record["justification_status"]),
                "justification_criterion": (
                    "none" if mapping_unresolved == "Y"
                    and record["justification_criterion"] == "unclear"
                    else record["justification_criterion"]),
                "source_cited": ("none" if mapping_unresolved == "Y"
                                 and record["source_cited"] == "unclear"
                                 else record["source_cited"]),
                "sensitivity_to_J": record["sensitivity_to_J"],
                "sensitivity_evidence": record["sensitivity_evidence"],
                "source_url": "",
                "source_location": record["source_location"],
                "source_version_checked": record["source_version"],
                "justification_given": GIVEN_RULE[justification_class],
                "justification_quote": record["justification_quote"],
                "secondary_repetition": record["other_counts_reported"],
                "coding_notes": record["coding_notes"],
                "J_is_outer_replication": record["J_is_outer_replication"],
                "confounded_with": record["confounded_with"],
                "justification_type": justification_class,
                "double_coded": "N",
                "coding_reviewer": "HGS",
            }
    return coded


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
            or i in FULLTEXT_ACCESSIBLE_ROWS
        )
    }
    if accessible_rows != set(MANUAL):
        raise ValueError(
            "open-access coding keys do not match the open-access stratum: "
            f"accessible={sorted(accessible_rows)}, manual={sorted(MANUAL)}"
        )

    local = load_local_pdf_coding(base_rows)
    clash = sorted(set(local) & set(MANUAL))
    if clash:
        raise ValueError(f"a corpus row is coded twice: {clash}")

    coded = {**MANUAL, **local}
    rows = []
    for corpus_row, base in enumerate(base_rows, 1):
        manual = coded.get(corpus_row, inaccessible_defaults(base))
        if corpus_row in coded:
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

    # Primary reporting rule first (PRIMARY_RULE): unclear counts as not justified.
    # The disaggregated four-way class follows immediately below, so a reader who
    # rejects the rule can reconstruct any other split from the same table.
    for category in ["justified", "not_justified"]:
        count = sum(row["justification_primary"] == category for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census, primary rule",
                "variable": "justification_primary",
                "category": category,
                "count": count,
                "denominator": len(eligible),
                "percentage": round(100 * count / len(eligible), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    # Sensitivity to J under the same rule: an unresolvable report is not evidence
    # that the paper varied J, so it is counted with the papers that did not.
    for category in ["yes", "no"]:
        if category == "yes":
            count = sum(row["sensitivity_to_J"] == "yes" for row in eligible)
        else:
            count = sum(row["sensitivity_to_J"] in {"no", "unclear"}
                        for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census, primary rule",
                "variable": "sensitivity_to_J_primary",
                "category": category,
                "count": count,
                "denominator": len(eligible),
                "percentage": round(100 * count / len(eligible), 4),
                "wilson_95_low": low,
                "wilson_95_high": high,
            }
        )

    # The row-level class is already three-way: ``unclear_report`` is folded into
    # ``unjustified`` in ``normalise_class``.  The ``mapping_unresolved`` block
    # below records how many of the unjustified records got there that way, so
    # the merge stays auditable from the summary alone.
    for category in ["explicit", "implicit_or_convention", "unjustified"]:
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

    # How many unjustified records reached that class through an unresolvable
    # count-to-source mapping rather than through a silent report.
    for category in ["Y", "N"]:
        count = sum(row["mapping_unresolved"] == category for row in eligible)
        low, high = wilson(count, len(eligible))
        summary.append(
            {
                "analysis_set": "accessible numeric-J census",
                "variable": "mapping_unresolved",
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
    print("primary rule (unclear counts as not justified):")
    for row in summary:
        if row["variable"] == "justification_primary":
            print(
                f"  {row['category']}: {row['count']}/{row['denominator']} "
                f"({row['percentage']:.2f}%, Wilson 95% CI "
                f"{row['wilson_95_low']:.2f}% to {row['wilson_95_high']:.2f}%)"
            )
    print("disaggregated class:")
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
