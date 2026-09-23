# Multiplicity: adjusting ten p-values by Holm, Hochberg and Benjamini-Hochberg.
#
# 1. DECIDE WHAT NEEDS ADJUSTING BEFORE LOOKING. One prespecified primary outcome needs none. A family of
#    secondary outcomes, subgroups or arms that will each support a claim does, and the family is named in
#    the protocol.
#
# 2. CHOOSE THE ERROR RATE, THEN THE METHOD.
#    - Family-wise error rate (the chance of ANY false claim): Holm, uniformly more powerful than
#      Bonferroni and valid under any dependence; Hochberg, more powerful again but valid only under
#      independence or positive dependence.
#    - False discovery rate (the expected share of false claims among the claims made): Benjamini-Hochberg,
#      for screening many outcomes where some false leads are acceptable.
#
# 3. REPORT ADJUSTED p-VALUES (p.adjust), and the unadjusted ones beside them, with the family and the
#    method named. Gatekeeping or graphical procedures (Bretz et al.) handle ordered hypotheses and are
#    a separate choice.

d <- read.csv("fixture.csv")
ys <- paste0("y", 1:10)
p <- sapply(ys, function(v) t.test(d[[v]][d$treated == 1], d[[v]][d$treated == 0])$p.value)
diffs <- sapply(ys, function(v) mean(d[[v]][d$treated == 1]) - mean(d[[v]][d$treated == 0]))
holm <- p.adjust(p, "holm")
hoch <- p.adjust(p, "hochberg")
bh <- p.adjust(p, "BH")

cat(sprintf("rejected at 0.05: unadjusted %d, Holm %d, Hochberg %d, Benjamini-Hochberg %d\n",
            sum(p < 0.05), sum(holm < 0.05), sum(hoch < 0.05), sum(bh < 0.05)))

cat("\n--- HARNESS ---\n")
for (j in 1:10) cat(sprintf("p_%d=%.10f\np_holm_%d=%.10f\np_hochberg_%d=%.10f\np_bh_%d=%.10f\n", j, p[j], j, holm[j], j, hoch[j], j, bh[j]))
cat(sprintf("difference_1=%.10f\n", diffs[1]))
cat(sprintf("rejected_unadjusted=%d\nrejected_holm=%d\nrejected_hochberg=%d\nrejected_bh=%d\n",
            sum(p < 0.05), sum(holm < 0.05), sum(hoch < 0.05), sum(bh < 0.05)))
