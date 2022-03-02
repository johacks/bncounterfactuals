set.seed(40)
bn <- NULL
if (!requireNamespace("bnlearn", quietly = TRUE)) {
  stop("bnlearn package needed to execute this functionality")
}
download.file("https://www.bnlearn.com/bnrepository/child/child.rda",
  "child.rda", "auto",
  quiet = TRUE
)
load("child.rda") # Load the popular CHILD Bayesian Network
file.remove("child.rda")
evidence <- data.frame(
  LVHReport = factor(x = "yes", levels = dimnames(bn$LVHreport$prob)[[1]]),
  LowerBodyO2 = factor(x = "5-12", levels = dimnames(bn$LowerBodyO2$prob)[[1]]),
  CO2Report = factor(x = "<7.5", levels = dimnames(bn$CO2Report$prob)[[1]]),
  XrayReport = factor(x = "Oligaemic", levels = dimnames(bn$XrayReport$prob)[[1]])
)
outcome <- predict(
  object = bn, node = bn$Disease$node, data = evidence, method = "bayes-lw"
)
print(paste("Outcome: ", outcome))
expected <- factor("TGA", levels = levels(outcome))
predict_f <- bnlearn_predict_wrapper(bn, bn$Disease$node, "bayes-lw", n = 5000)
r <- bfs_cfx(predict_f, evidence, outcome, expected)
print(data.table::rbindlist(r, fill = T))
