# R oracle for the TISCA v2 Python <-> CRAN MCS parity test.
#
# Runs the *actual* CRAN 'MCS' package on a loss matrix supplied as CSV and emits:
#   1. per-model MCS results (avg loss, H0 p-value, MCS p-value, included/excluded)
#   2. the exact bootstrap resample indices CRAN uses internally, so the Python
#      implementation can be driven with identical resamples (bit-exact parity).
#
# Usage: Rscript mcs_parity_oracle.R <loss.csv> <out_prefix> <alpha> <B> <statistic> <k> <seed>
# The loss CSV has rows = replications, columns = models, NO row/col header expected,
# numeric only (a single header row is stripped if present).
suppressMessages(requireNamespace("MCS"))

args <- commandArgs(trailingOnly = TRUE)
loss_file  <- args[1]
out_prefix <- args[2]
alpha      <- as.numeric(args[3])
B          <- as.integer(args[4])
statistic  <- args[5]
k          <- as.integer(args[6])
seed       <- as.integer(args[7])

mL <- as.matrix(read.csv(loss_file, header = FALSE))
Trows <- nrow(mL)

# --- 1. The exact bootstrap indices CRAN's MCSprocedure uses. ----------------
# MCSprocedure does: set.seed(seed); mIndices = GetIndices(iT, k, B).
set.seed(seed)
mIndices <- MCS:::GetIndices(Trows, k, B)
# Emit 0-based indices so Python can index numpy arrays directly.
write.table(mIndices - 1, file = paste0(out_prefix, "_indices.csv"),
            row.names = FALSE, col.names = FALSE, sep = ",")

# --- 2. The CRAN result. ------------------------------------------------------
set.seed(seed)
res <- MCS::MCSprocedure(Loss = mL, alpha = alpha, B = B, statistic = statistic,
                         k = k, verbose = FALSE, seed = seed)

tab <- res@show  # matrix with columns Avg.Loss, p-Value for H0, MCS p-Value
out <- data.frame(
    model    = rownames(tab),
    avg_loss = tab[, "Avg.Loss"],
    p_H0     = tab[, "p-Value for H_{0,M_k}"],
    p_mcs    = tab[, "MCS p-Value"]
)
write.csv(out, file = paste0(out_prefix, "_result.csv"), row.names = FALSE)

# Included / excluded model lists.
write.csv(data.frame(included = res@Info$included),  file = paste0(out_prefix, "_included.csv"),  row.names = FALSE)
write.csv(data.frame(excluded = res@Info$excluded),  file = paste0(out_prefix, "_excluded.csv"),  row.names = FALSE)
cat("DONE\n")
