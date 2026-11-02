#!/usr/bin/env Rscript
# Readout for prereg X: planted-truth gates (power at OR 1.5 and 2, false positives, confound arms), then, if the
# real run exists, the frozen label. Rules are prereg X's, applied mechanically.
# Usage, from analysis/: Rscript breeder_diffusion_readout.R
# Output: ../data/breeder_diffusion_readout.json
suppressPackageStartupMessages(library(jsonlite))
rd <- function(arm) {
  m <- do.call(rbind, lapply(Sys.glob(sprintf("../data/breeder_diffusion_sim_%s_*.json", arm)), function(f) fromJSON(f)$reps))
  colnames(m) <- c("beta", "se", "or", "lo", "hi", "p", "n_events", "n_rows"); stopifnot(nrow(m) == 100); m
}
flag <- function(m) sum(m[, "p"] < 0.05 & m[, "beta"] > 0)
g <- sapply(c("or1", "or1.5", "or2", "conf0.5", "conf1"), function(a) {
  m <- rd(a); c(flagged = flag(m), any_sig = sum(m[, "p"] < 0.05), median_or = median(m[, "or"]), median_events = median(m[, "n_events"]))
})
print(round(t(g), 3))
underpowered <- g["flagged", "or1.5"] < 80
invalid <- g["any_sig", "or1"] > 8
confounded <- g["flagged", "conf0.5"] > 8
cat(sprintf("gates: underpowered at OR 1.5 %s (power at OR 2: %d/100); invalid %s; confound arm fools the test %s (%d/100 at gamma 0.5, %d/100 at 1)\n",
            underpowered, g["flagged", "or2"], invalid, confounded, g["flagged", "conf0.5"], g["flagged", "conf1"]))
out <- list(gates = g, underpowered = underpowered, invalid = invalid, confounded = confounded)

real_f <- "../data/breeder_diffusion_real.json"
if (file.exists(real_f)) {
  r <- fromJSON(real_f); m <- r$main
  lab <- if (invalid) "invalid" else if (m$p < 0.05 && m$or > 1) {
    if (confounded) "adoption clustered on the breeder network; copying and shared background cannot be told apart"
    else "followed the network (consistent with copying)"
  } else if (m$p >= 0.05 && m$hi < 1.5 && !underpowered) "no network effect detected" else "unresolved"
  cat(sprintf("real: OR %.2f (%.2f-%.2f), p %.4f, %d adoptions in %d breeder-years; label: %s\n", m$or, m$lo, m$hi, m$p,
              m$n_events, m$n_rows, lab))
  out$real <- list(main = m, label = lab, registrant = r$registrant, walkeriana = r$walkeriana, binary = r$binary)
}
write_json(out, "../data/breeder_diffusion_readout.json", digits = 6, auto_unbox = TRUE)
