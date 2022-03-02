# Get counterfactuals / sufficient explanations with Bayesian Networks

################### Definitions
`%>%` <- magrittr::`%>%`

LABEL_TRUE <- 0L
LABEL_EXP <- 1L
LABEL_OTHER <- 2L
LABEL_NULL <- 3L

################### Utility functions

# Get all possible parents of a subset of variables
enumerate_parents <- coro::generator(function(label_subset_mask) {
  for (i in 1:length(label_subset_mask)) {
    if (!label_subset_mask[i]) {
      label_subset_mask[i] <- TRUE
      coro::yield(label_subset_mask)
      label_subset_mask[i] <- FALSE
    }
  }
})

# Get all possible children of a subset of variables
enumerate_children <- coro::generator(function(label_subset_mask) {
  for (i in 1:length(label_subset_mask)) {
    if (label_subset_mask[i]) {
      label_subset_mask[i] <- FALSE
      coro::yield(label_subset_mask)
      label_subset_mask[i] <- TRUE
    }
  }
})

possible_cfxs <- function(label_dict, bool_mask.s, bool_mask.c, evidence,
                          hypothetical_values) {
  # Values of original evidence for current subset of variables
  values.s <- evidence[1, bool_mask.s]
  # Possible different values of rest of variables
  values.c <- expand.grid(... = hypothetical_values[bool_mask.c])
  # Combine both above
  p_cfxs <- evidence[FALSE, ] # Same columns, but empty
  p_cfxs[1:nrow(values.c), bool_mask.c] <- values.c
  p_cfxs[, bool_mask.s] <- values.s
  return(p_cfxs)
}

# For a subset of variables, compute modes and associated label
compute_modes_and_label <- function(label_dict, label_subset_mask, evidence,
                                    outcome, expected, predict_f, p_cfxs) {
  # Get prediction for hypothetical situations
  prediction <- unlist(predict_f(p_cfxs))
  # Annotate label
  if (all(prediction == outcome)) {
    label_dict[[label_subset_mask]] <- LABEL_TRUE
  } else if (all(prediction == expected)) {
    label_dict[[label_subset_mask]] <- LABEL_EXP
  } else {
    label_dict[[label_subset_mask]] <- LABEL_OTHER
  }
  return(prediction)
}

# All possible different values that each variable in the evidence can have,
# except for the current one
generate_hypothetical_values <- function(evidence) {
  evidence %>% Map(f = function(val) {
    Filter(x = levels(val), f = function(val_level) val_level != val)
  })
}

# All elements of l1 are in l2
is_subset <- function(l1, l2) {
  all(purrr::map2_lgl(
    names(l1), l1, function(name, e1) !is.null(e2 <- l2[[name]]) && e1 == e2
  ))
}

# Create a label dictionary for use in bfs_sfx/cfx algorithms
create_label_dict <- function() {
  r2r::hashmap(default = LABEL_NULL)
}

################### Public functions

bfs_sfx_cfx_aux <- function(predict_f, evidence, outcome, expected, sfx_only = F) {
  # Start point is all variables included, annotate the label
  label_subset_mask <- bit::as.bit(rep(TRUE, length(evidence)))
  label_dict <- create_label_dict()
  # Possible different values of each variable in the evidence
  hypothetical_values <- generate_hypothetical_values(evidence)
  # Queue children
  pot_s <- list()
  queue_sfx <- dequer::as.queue(label_subset_mask)
  if (!sfx_only) {
    queue_cfx <- dequer::as.deque(label_subset_mask)
    cfxs <- list()
  }
  while (length(queue_sfx) != 0) {
    label_subset_mask <- dequer::pop(queue_sfx)
    all_parents_true <- TRUE
    coro::loop(for (p in enumerate_parents(label_subset_mask)) {
      if (!(all_parents_true <- label_dict[[p]] == LABEL_TRUE)) break
    })
    if (label_dict[[label_subset_mask]] == LABEL_NULL && all_parents_true) {
      suff_mask_bool <- as.logical(label_subset_mask)
      cfx_mask_bool <- !suff_mask_bool
      pot_cfx <- possible_cfxs(
        label_dict, suff_mask_bool,
        cfx_mask_bool, evidence, hypothetical_values
      )
      prediction <- compute_modes_and_label(
        label_dict, label_subset_mask, evidence, outcome,
        expected, predict_f, pot_cfx
      )
      if (label_dict[[label_subset_mask]] == LABEL_TRUE) {
        pot_s[[length(pot_s) + 1]] <- label_subset_mask
        coro::loop(for (c in enumerate_children(label_subset_mask)) {
          dequer::pushback(queue_sfx, c)
        })
      } else if (!sfx_only) {
        rows_exp <- prediction == expected
        if (nrow(match <- pot_cfx[rows_exp, cfx_mask_bool, drop = F]) > 0) {
          for (i in 1:nrow(match)) {
            cfxs[[(length(cfxs) + 1)]] <- as.list(match[i, , drop = F])
          }
        }
        if (label_dict[[label_subset_mask]] != LABEL_EXP) {
          coro::loop(for (c in enumerate_children(label_subset_mask)) {
            dequer::pushback(queue_cfx, c)
          })
        }
      }
    }
  }
  # Get resulting sufficient explanations
  sfxs <- list()
  for (label_subset_mask in pot_s) {
    all_children_not_true <- TRUE
    coro::loop(for (c in enumerate_children(label_subset_mask)) {
      if (!(label_dict[[c]] != LABEL_TRUE -> all_children_not_true)) break
    })
    if (all_children_not_true) {
      idx <- as.logical(label_subset_mask)
      sfxs[[length(sfxs) + 1]] <- as.list(evidence[1, idx, drop = FALSE])
    }
  }
  if (!sfx_only) {
    cfxs <- get_remaining_cfx(
      predict_f, evidence, outcome, expected,
      label_dict, queue_cfx, cfxs
    )
    return(list(sfxs = sfxs, cfxs = cfxs))
  } else {
    return(sfxs)
  }
}

#' Get sufficient and counterfactual explanations for a prediction
#'
#' @param predict_f A function that takes in a dataframe of evidence and outputs
#'                  a list with the prediction for each row.
#' @param evidence A dataframe of one row where all columns are of type factor.
#' @param outcome The resulting prediction for the evidence.
#' @param expected The expected prediction of the evidence (only for
#'                 counterfactuals).
#'
#' @return Sufficient explanations as a list of named lists.
#' @export
#'
#' @example inst/examples/bfs_sfx_cfx.R
bfs_sfx_cfx <- function(predict_f, evidence, outcome, expected) {
  bfs_sfx_cfx_aux(predict_f, evidence, outcome, expected, sfx_only = F)
}

#' Get sufficient explanations for a prediction
#'
#' @param predict_f A function that takes in a dataframe of evidence and outputs
#'                  a list with the prediction for each row.
#' @param evidence A dataframe of one row where all columns are of type factor.
#' @param outcome The resulting prediction for the evidence.
#' @param expected The expected prediction of the evidence
#'
#' @return Sufficient explanations as a list of named lists.
#' @export
#'
#' @example inst/examples/bfs_sfx.R
bfs_sfx <- function(predict_f, evidence, outcome, expected) {
  bfs_sfx_cfx_aux(predict_f, evidence, outcome, expected, sfx_only = T)
}

#' Get counterfactual explanations for a prediction
#'
#' @param predict_f A function that takes in a dataframe of evidence and outputs
#'                  a list with the prediction for each row.
#' @param evidence A dataframe of one row where all columns are of type factor.
#' @param outcome The resulting prediction for the evidence.
#' @param expected The expected prediction of the evidence
#'
#' @return Counterfactual explanations as a list of named lists.
#' @export
#'
#' @example inst/examples/bfs_cfx.R
bfs_cfx <- function(predict_f, evidence, outcome, expected) {
  # Start point is all variables included, annotate the label
  label_subset_mask <- bit::as.bit(rep(TRUE, length(evidence)))
  label_dict <- create_label_dict()
  cfxs <- list() # Counterfactual explanation set initially empty
  queue <- dequer::as.queue(label_subset_mask)
  get_remaining_cfx(
    predict_f, evidence, outcome, expected, label_dict, queue, cfxs
  )
}

# Queue processing for counterfactual explanations
get_remaining_cfx <- function(predict_f, evidence, outcome, expected,
                              label_dict, queue, cfxs) {
  # Possible different values of each variable in the evidence
  hypothetical_values <- generate_hypothetical_values(evidence)

  while (length(queue) != 0) {
    label_subset_mask <- dequer::pop(queue)
    no_parent_exp <- TRUE
    coro::loop(for (p in enumerate_parents(label_subset_mask)) {
      if (!(no_parent_exp <- label_dict[[p]] != LABEL_EXP)) break
    })
    if (label_dict[[label_subset_mask]] == LABEL_NULL && no_parent_exp) {
      suff_mask_bool <- as.logical(label_subset_mask)
      cfx_mask_bool <- !suff_mask_bool
      pot_cfx <- possible_cfxs(
        label_dict, suff_mask_bool, cfx_mask_bool,
        evidence, hypothetical_values
      )
      # No already seen cfx can be a subset of a new one
      rows <- purrr::map_lgl(1:nrow(pot_cfx), function(i) {
        !any(purrr::map_lgl(cfxs, function(l) is_subset(l, pot_cfx[i, ])))
      })
      pot_cfx <- pot_cfx[rows, ]
      if (nrow(pot_cfx) > 0) {
        # Compute modes and label and keep track of successful counterfactuals
        prediction <- compute_modes_and_label(
          label_dict, label_subset_mask, evidence, outcome, expected, predict_f,
          pot_cfx
        )
        rows_exp <- prediction == expected
        if (nrow(match <- pot_cfx[rows_exp, cfx_mask_bool, drop = F]) > 0) {
          for (i in 1:nrow(match)) {
            cfxs[[(length(cfxs) + 1)]] <- as.list(match[i, , drop = F])
          }
        }
        if (label_dict[[label_subset_mask]] != LABEL_EXP) {
          coro::loop(for (c in enumerate_children(label_subset_mask)) {
            dequer::pushback(queue, c)
          })
        }
      }
    }
  }
  return(cfxs)
}

#' Create a predict function for use of a bnlearn Bayesian network in the
#' search algorithms.
#'
#' @param bn Fitted Bayesian network
#' @param node Name of target node to be predicted
#' @param method Same options as in parameter 'method' in predict function of
#' bnlearn
#' @param threshold Minimum value of posterior probability of the MLE to be
#' the output of the algorithm. If parameter 'default' is null, the MLE will
#' be the output even if it doesn't surpass the threshold.
#' @param default Default predicted target value when threshold isn't exceeded
#' by the MLE
#' @param n Same options as in parameter 'n' in predict function of bnlearn
#'
#' @return Function to be used in algorithms such as bfs_sfx
#' @export
#'
#' @example inst/examples/bnlearn_predict_wrapper.R
bnlearn_predict_wrapper <- function(bn, node, method, threshold = NULL,
                                    default = NULL, n = 500) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("bnlearn package needed to execute this functionality")
  }
  levs <- dimnames(bn[[node]]$prob)[[1]]
  fact <- factor(levs, levels = levs)

  return(function(df) {
    if (method == "bayes-lw") {
      pred <- predict(
        object = bn, node = node,
        data = df, method = method, n = n, prob = TRUE
      )
    } else {
      pred <- predict(
        object = bn, node = node,
        data = df, method = method, prob = TRUE
      )
    }

    probs <- attr(pred, "prob")
    Map(function(i) {
      mle_i <- which.max(probs[, i])
      mle <- probs[mle_i, i]
      if (!is.null(threshold)) {
        if (mle > threshold) {
          return(fact[[mle_i]])
        } else if (is.null(default)) {
          return(fact[[mle_i]])
        } else {
          return(default)
        }
      }
      return(fact[[mle_i]])
    }, 1:length(pred))
  })
}


################### TESTS: pending revision

bfs_sfx_test <- function() {
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
  expected <- factor("PAIVS", levels = levels(outcome))
  predict_f <- bnlearn_predict_wrapper(bn, bn$Disease$node, "bayes-lw",
    threshold = 0.20, default = "TGA", n = 500000
  )
  r <- bfs_sfx(predict_f, evidence, outcome, expected)
  print(data.table::rbindlist(r, fill = T))
}

bfs_cfx_test <- function() {
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
  predict_f <- bnlearn_predict_wrapper(bn, bn$Disease$node, "bayes-lw", n = 500000)
  r <- bfs_cfx(predict_f, evidence, outcome, expected)
  print(data.table::rbindlist(r, fill = T))
}

bfs_sfx_cfx_test <- function() {
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
  predict_f <- bnlearn_predict_wrapper(bn, bn$Disease$node, "bayes-lw", n = 500000)
  r <- bfs_sfx_cfx(predict_f, evidence, outcome, expected)
  print("Sufficient explanations:")
  print(data.table::rbindlist(r[["sfxs"]], fill = T))
  print("Counterfactual explanations")
  print(data.table::rbindlist(r[["cfxs"]], fill = T))
}
