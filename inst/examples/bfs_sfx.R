set.seed(40)
download.file("https://www.bnlearn.com/bnrepository/child/child.rda",
  "child.rda", "auto",
  quiet = TRUE
)
load("child.rda") # Load CHILD Bayesian Network
evidence <- data.frame(
  LVHReport = factor(x = "yes", levels = dimnames(bn$LVHreport$prob)[[1]]),
  LowerBodyO2 = factor(x = "5-12", levels = dimnames(bn$LowerBodyO2$prob)[[1]]),
  CO2Report = factor(x = "<7.5", levels = dimnames(bn$CO2Report$prob)[[1]]),
  XrayReport = factor(x = "Oligaemic", levels = dimnames(bn$XrayReport$prob)[[1]])
)
outcome <- predict(
  object = bn, node = bn$Disease$node, data = evidence,
  method = "bayes-lw"
)
print(paste("Outcome: ", outcome))
expected <- factor("TGA", levels = levels(outcome))
predict_f <- function(df) {
  predict(object = bn, node = bn$Disease$node, data = df, method = "bayes-lw")
}
print(bfs_sfx(predict_f, evidence, outcome, expected))
