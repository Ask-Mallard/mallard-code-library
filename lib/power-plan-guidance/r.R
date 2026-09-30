# Structured planning guidance when an ordinary power calculation is invalid.

planning_guidance <- function(approach, why, required_inputs, next_step) {
  allowed <- c("precision", "fixed-design", "rule-of-thumb", "bespoke-simulation", "not-applicable")
  stopifnot(approach %in% allowed, nzchar(why), nzchar(required_inputs), nzchar(next_step))
  list(
    approach=approach,
    powerCalculated=FALSE,
    sampleSizeCalculated=FALSE,
    why=why,
    requiredInputs=strsplit(required_inputs, ";", fixed=TRUE)[[1]],
    nextStep=next_step
  )
}

d <- read.csv("fixture.csv", stringsAsFactors=FALSE)
plans <- lapply(seq_len(nrow(d)), function(i) planning_guidance(d$approach[i],d$why[i],d$required_inputs[i],d$next_step[i]))
for (plan in plans) {
  cat(sprintf("%s\n  Why: %s\n  Required: %s\n  Next: %s\n",plan$approach,plan$why,paste(plan$requiredInputs,collapse=", "),plan$nextStep))
}
cat("\n--- HARNESS ---\n")
cat(sprintf("route_count=%d\n",length(plans)))
cat(sprintf("numeric_sample_sizes=%d\n",sum(vapply(plans,function(x)x$sampleSizeCalculated,logical(1)))))
cat(sprintf("routes_with_next_step=%d\n",sum(vapply(plans,function(x)nzchar(x$nextStep),logical(1)))))
