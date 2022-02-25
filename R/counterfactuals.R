# Get counterfactuals / sufficient explanations with Bayesian Networks

# Libraries
`%>%` <- magrittr::`%>%`

LABEL_TRUE <- 0L
LABEL_EXP <- 1L
LABEL_OTHER <- 2L
LABEL_NULL <- 3L

# Get all possible parents of a subset of variables
enumerate_parents <- coro::generator(function (label_subset_mask) {
  for (i in 1:length(label_subset_mask))
    if (!label_subset_mask[i]){
      label_subset_mask[i] = TRUE
      coro::yield(label_subset_mask)
      label_subset_mask[i] = FALSE
    }
})

# Get all possible children of a subset of variables
enumerate_children <- coro::generator(function (label_subset_mask) {
  for (i in 1:length(label_subset_mask)){
    if (label_subset_mask[i]){
      label_subset_mask[i] = FALSE
      coro::yield(label_subset_mask)
      label_subset_mask[i] = TRUE
    }
  }

})

# For a subset of variables, compute modes and associated label
compute_modes_and_label <- function(label_dict, label_subset_mask, evidence,
                                    target, expected, predict_f,
                                    hypothetical_values) {

  bool_mask.suff <- as.logical(label_subset_mask)
  bool_mask.cntrfact <- !bool_mask.suff
  # Values of original evidence for current subset of variables
  values.suff <- evidence[1, bool_mask.suff]
  # Possible different values of rest of variables
  values.cntrfact <- expand.grid(...=hypothetical_values[bool_mask.cntrfact])
  # Combine both above
  values.mix <- evidence[FALSE, ] # Same columns, but empty
  values.mix[1:nrow(values.cntrfact), bool_mask.cntrfact] <- values.cntrfact
  values.mix[, bool_mask.suff] <- values.suff
  # Get prediction for hypothetical situations
  prediction <- unlist(predict_f(values.mix))
  # Annotate label
  if (all(prediction == target))
    label_dict[[label_subset_mask]] = LABEL_TRUE
  else if (all(prediction == expected))
    label_dict[[label_subset_mask]] = LABEL_EXP
  else
    label_dict[[label_subset_mask]] = LABEL_OTHER
}

# All possible different values that each variable in the evidence can have,
# except for the current one
generate_hypothetical_values <- function(evidence) {
  evidence %>% Map(f = function(val) {
      Filter(x = levels(val), f = function(val_level) val_level != val)
    })
}

#' Get sufficient explanations for a prediction
#'
#' @param predict_f A function that takes in a dataframe of evidence and outputs
#'                  a list with the prediction for each row.
#' @param evidence A dataframe of one row where all columns are of type factor.
#' @param target The resulting prediction for the evidence.
#' @param expected The expected prediction of the evidence (only for
#'                 counterfactuals).
#'
#' @return List where each element is a sufficient explanation of the prediction.
#' @export
#'
#' @examples
#' library(bnlearn)
#' set.seed(40)
#' download.file("https://www.bnlearn.com/bnrepository/child/child.rda",
#'               "child.rda", "auto", quiet = TRUE)
#' load("child.rda")  # Load CHILD Bayesian Network
#' evidence = data.frame(
#'   LVHReport=factor(x = "yes", levels = dimnames(bn$LVHreport$prob)[[1]]),
#'   LowerBodyO2=factor(x = "5-12", levels = dimnames(bn$LowerBodyO2$prob)[[1]]),
#'   CO2Report=factor(x = "<7.5", levels = dimnames(bn$CO2Report$prob)[[1]]),
#'   XrayReport=factor(x = "Oligaemic", levels = dimnames(bn$XrayReport$prob)[[1]])
#' )
#' target <- predict(object = bn, node = bn$Disease$node, data = evidence,
#'                  method = "bayes-lw")
#' print(target)
#' expected <- factor("PAIVS", levels = levels(target))
#' predict_f <- function(df) {
#'   predict(object = bn, node = bn$Disease$node, data = df, method = "bayes-lw")}
#' print(bfs_sfx(predict_f, evidence, target, expected))
bfs_sfx <- function(predict_f, evidence, target, expected) {
  # Start point is all variables included, annotate the label
  label_subset_mask = bit::as.bit(rep(TRUE, length(evidence)))
  label_dict <- r2r::hashmap(default = LABEL_NULL)
  label_dict[[label_subset_mask]] = LABEL_TRUE
  # Possible different values of each variable in the evidence
  hypothetical_values <- generate_hypothetical_values(evidence)
  # Queue children
  pot_s <- list(label_subset_mask)
  pending <- dequer::queue()
  coro::loop(for (c in enumerate_children(label_subset_mask))
       dequer::pushback(pending, c))
  # Process pending children until finished
  while (length(pending) != 0) {
    label_subset_mask <- dequer::pop(pending)
    all_parents_true <- TRUE
    coro::loop(for (p in enumerate_parents(label_subset_mask)) {
      all_parents_true <- all_parents_true && label_dict[[p]] == LABEL_TRUE
      if (!all_parents_true) break
    })
    if (label_dict[[label_subset_mask]] == LABEL_NULL && all_parents_true) {
      compute_modes_and_label(label_dict, label_subset_mask, evidence, target,
                              expected, predict_f, hypothetical_values)
      if (label_dict[[label_subset_mask]] == LABEL_TRUE) {
        pot_s[[length(pot_s)+1]] <- label_subset_mask
        coro::loop(for (c in enumerate_children(label_subset_mask))
                   dequer::pushback(pending, c))
      }
    }
  }
  # Get resulting sufficient explanations
  sufficient_exps <- list()
  for (label_subset_mask in pot_s) {
    all_children_not_true = TRUE
    coro::loop(for (c in enumerate_children(label_subset_mask)) {
      all_children_not_true <- all_children_not_true && label_dict[[c]] != LABEL_TRUE
      if (!all_children_not_true) break
    })
    if (all_children_not_true) {
      idx = as.logical(label_subset_mask)
      sufficient_exps[[length(sufficient_exps)+1]] <- names(evidence)[idx]
    }
  }
  return(sufficient_exps)
}

bfs_sfx_test <- function() {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    warning("bnlearn package needed to execute this test")
    return()
  }
  bn.predict <- utils::getFromNamespace("predict", "bnlearn")
  set.seed(40)
  bn <- NULL
  load("child.rda")  # Load the popular CHILD Bayesian Network
  evidence = data.frame(
    LVHReport=factor(x = "yes", levels = dimnames(bn$LVHreport$prob)[[1]]),
    LowerBodyO2=factor(x = "5-12", levels = dimnames(bn$LowerBodyO2$prob)[[1]]),
    CO2Report=factor(x = "<7.5", levels = dimnames(bn$CO2Report$prob)[[1]]),
    XrayReport=factor(x = "Oligaemic", levels = dimnames(bn$XrayReport$prob)[[1]]))
  target <- bn.predict(
    object = bn, node = bn$Disease$node, data = evidence,method = "bayes-lw")
  expected <- factor("PAIVS", levels = levels(target))
  predict_f <- function(df) {
    bn.predict(
      object = bn, node = bn$Disease$node, data = df, method = "bayes-lw")}

  print(bfs_sfx(predict_f, evidence, target, expected))
}

#' Title
#'
#' @param bn
#' @param node
#' @param method
#' @param threshold
#' @param default
#' @param n
#'
#' @return
#' @export
#'
#' @examples
bnlearn_predict_wrapper <- function(bn, node, method, threshold = NULL,
                                    default = NULL, n = 500){
  return (function(df) {
    levs <- dimnames(red.k2.fit[[node]]$prob)[[1]]
    fact <- factor(levs, levels = levs)
    if (method == "bayes-lw") {
      pred <- predict(object = red.k2.fit, node = red.k2.fit$Muerto$node,
                      data = df, method = method, n=n, prob = TRUE)
    }
    else pred <- predict(object = red.k2.fit, node = red.k2.fit$Muerto$node,
                         data = df, method = method, prob = TRUE)

    probs <- attr(pred, "prob")
    Map(function(i) {
      mle_i <- which.max(probs[,i])
      mle <- probs[max, i]
      if (!is.null(threshold)) {
        if (mle > threshold) return(fact[[mle_i]])
        else if (is.null(default)) return(fact[[mle_i]])
        else return(default)
      }
      return(fact[[mle_i]])
    }, 1:length(pred))
  })
}

debug_p <- function(){
  library(bnlearn)
  load("../Notebooks/datos-estaticos.rdata")
  load("../Notebooks/red-estatica.rdata")
  comorbilidades <-
    list("Fumador",
         "Cardio",
         "Pulmonar",
         "Diabetes",
         "Renal",
         "Neuro",
         "Onco",
         "HTA")
  evidencia <- datos.red[FALSE, ]
  evidencia <- datos.red[1, unlist(comorbilidades)]
  resultado <- predict(object = red.k2.fit, node = red.k2.fit$Muerto$node,
                       data = evidencia, method = "bayes-lw", n=1000)
  predict_f <- bnlearn_predict_wrapper(red.k2.fit, red.k2.fit$Muerto$node,
                                       "bayes-lw", 0.7)
  esperado = "S"
  res <- bfs_sfx(predict_f = predict_f, evidence = evidencia,
                 target = resultado, expected = esperado)
  print("Evidencia: ")
  print(evidencia)
  print(res)
}

