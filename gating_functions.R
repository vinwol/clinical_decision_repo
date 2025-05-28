pbetaMix <- function(x, par, weights, lower.tail=TRUE) {
  ret <- sum(weights * pbeta(x, par[, 1], par[, 2], lower.tail=lower.tail))
  return(ret)
}

pbetaMix <- Vectorize(pbetaMix, vectorize.args="x")

getBetamixPost <- function(x, n, par, weights) {
  # Check the format.
  stopifnot(is.matrix(par),
            is.numeric(par),
            identical(ncol(par), 2L),
            all(par > 0),
            identical(nrow(par), length(weights)),
            all(weights > 0))
  
  # Renormalize weights.
  weights <- weights / sum(weights)
  
  # Compute updated parameters.
  postPar <- par
  postPar[, 1] <- postPar[, 1] + x
  postPar[, 2] <- postPar[, 2] + n - x
  postParProb <- postPar[, 1] / (postPar[, 1] + postPar[, 2])
  
  # Compute updated mixture probabilities.
  tmp <- exp(stats::dbinom(x, size = n, prob = postParProb, log = TRUE) +
               stats::dbeta(postParProb, par[, 1], par[, 2], log = TRUE) -
               stats::dbeta(postParProb, postPar[, 1], postPar[, 2], log = TRUE))
  
  postWeights <- weights * tmp / sum(weights * tmp)
  
  return(list(par=postPar, weights=postWeights))
}

postprob <- function(x, n, p, parE=c(1, 1), weights, betamixPost, log.p=FALSE) {
  
  if(missing(betamixPost)) {
    # if parE is a vector => situation where there is only one component.
    if(is.vector(parE)) {
      # Check that it has exactly two entries.
      stopifnot(identical(length(parE), 2L))
      # Transpose to matrix with one row.
      parE <- t(parE)
    }
    
    # If prior weights of the beta mixture are not supplied.
    if(missing(weights)) {
      weights <- rep(1, nrow(parE))
      # (don't need to be normalized, this is done in getBetamixPost)
    }
    
    # Compute updated parameters.
    betamixPost <- getBetamixPost(x=x, 
                                  n=n,
                                  par=parE,
                                  weights=weights)
  }
  
  # Compute the survival function at p, i.e. 1 - cdf at p:
  ret <- with(betamixPost,
              pbetaMix(x=p, par=par, weights=weights, lower.tail=FALSE))
  
  if(log.p) {
    return(log(ret))
  } else {
    return(ret)
  }
}

postprob <- Vectorize(postprob, vectorize.args="x")

is_not_valid_number <- function(var) {
  !is.finite(var) || is.na(var) || is.nan(var)
}

percent_conditional_round_vectorized <- function(vec) {
  vec <- vec * 100
  new_vec <- c()
  for (elem in vec) {
    if (elem < 1) {
      new_elem <- round(elem,2)
      new_vec <- append(new_vec,new_elem)
    } else {
      new_elem <- round(elem,1)
      new_vec <- append(new_vec,new_elem)
    }
  }  
  return(new_vec)
}

create_futility_data <- function(max_responders,sample_size) {
  x <- seq(0,1,0.01) 
  binomial_prob <- stats::pbinom(max_responders, 
                                 sample_size, 
                                 x, 
                                 lower.tail = FALSE, 
                                 log.p = FALSE)
  r <- rep(max_responders,length(binomial_prob))
  n <- rep(sample_size,length(binomial_prob))
  orr <- r/n
  df <- data.frame(x=x,
                   y=(1-binomial_prob),
                   r=r,
                   n=n,
                   observed_orr = orr)
  return(df)
}


# Operational characteristics for decision gating.
# The probability of Go or Stop follows a Binomial distribution.
# binomial(successes = number of responders, size = number of trials or patients, 
# prob = True RR)
OCprob <- function(x, 
                   N, 
                   resp.go, 
                   resp.fut, 
                   prob) {
  p.go <- round(1 - stats::pbinom(q=resp.go-1, size=N, prob=prob), 4)
  # resp.fut = number of sucesses in futile context, which shares same prob of 
  # success in " Go" case
  p.stop <- round(stats::pbinom(q=resp.fut, size=N, prob=prob), 4) 
  df <- data.frame(prob=prob, p.go=p.go, p.stop = p.stop)
  return(df)
}

get_responder_go <- function(n, 
                             orr_tpp_go,
                             confidence_go) {
  Nresp <- 0:n
  pp_go <- postprob(Nresp, # number of successes
                    n, # number of patients
                    orr_tpp_go, # threshold
                    parE = c(1, 1), # beta parameters, here uniform prior
                    log.p = FALSE) # log of the probability? (default: FALSE)
  min_responders_go <- min(Nresp[pp_go >= confidence_go])
  return(min_responders_go)
}

get_responder_stop <- function(n, 
                               orr_tpp_stop,
                               confidence_stop) {
  Nresp <- 0:n
  pp_stop <- postprob(Nresp, 
                      n, 
                      orr_tpp_stop, 
                      parE = c(1, 1), 
                      log.p = FALSE)
  max_responders_stop <- max(Nresp[pp_stop < confidence_stop])
  return(max_responders_stop)
}
# Dataframe to generate posterior for a number of responders (e.g. 0:12), given an
# uninformative prior Beta(alpha = 1, beta = 1)
# 1- pbeta(True PP, 0:12 + alpha, beta + n.pat - 0:12)
# We take from the probs table that the min True RR that satifies a "Go!" 
# corresponds to this number of events or responders.
# We take from the probs table that the min True RR that satifies a "Stop!" 
# corresponds to this number of events or responders, could also be satisfied with 
# probs[(probs$prob)<=conf.stop,][1,1]
generate_posterior_probability_data <- function(num_patients, min_posterior_prob) {
  
  responders <- seq(0, num_patients)
  
  # compute the cumulative distribution function (CDF) of the beta distribution.
  # q: the quantile at which to compute the CDF.
  # The quantile is given by MinPP (0.45).
  # shape1: the first shape parameter of the beta distribution (usually denoted as alpha).
  # The first shape parameter is seq(0, n.pat) + 1, which creates a sequence from 1 to n.pat + 1.
  # shape2: the second shape parameter of the beta distribution (usually denoted as beta).
  # The second shape parameter is n.pat - seq(0, n.pat) + 1, which creates a sequence from n.pat + 1 to 1.
  # ncp: (optional) the non-centrality parameter (not used in this case).
  # This means that for each value in the sequence of responders from 0 to n.pat (40 in this case), 
  # you are computing the CDF at MinPP for a beta distribution parameterized by 
  # corresponding shape1 and shape2 values.
  # Thus, the pbeta function evaluates the beta CDF at these parameter combinations, 
  # and then you subtract the result from 1 to get the complementary probability.
  # We subtract the cumulative distribution function (CDF) of the beta distribution 
  # from one to obtain the complementary probability, which represents the upper tail 
  # probability of the corresponding beta distribution. 
  # This is often done in hypothesis testing and other statistical procedures to 
  # determine the probability of observing a more extreme value than the given quantile.
  posterior_prob <- 1 - stats::pbeta(min_posterior_prob, 
                                     responders + 1, 
                                     num_patients - responders + 1)
  
  # Create dataframe with responders and their corresponding probabilities.
  posterior_probability_data <- data.frame(responders = responders, 
                                           posterior_prob = posterior_prob)
  
  return(posterior_probability_data)
}

# Operational characteristics for decision gating.
# The probability of Go or Stop follows a Binomial distribution.
# binomial(successes = number of responders, size = number of trials or patients, prob = True RR)
oc_prob <- function(x, 
                    N, 
                    resp.go, 
                    resp.fut, 
                    prob) {
  p.go <- round(1 - stats::pbinom(q=resp.go-1, size=N, prob=prob), 4)
  # resp.fut = number of sucesses in futile context, which shares same prob of 
  # success in " Go" case
  p.stop <- round(stats::pbinom(q=resp.fut, size=N, prob=prob), 4) 
  df <- data.frame(prob=prob, p.go=p.go, p.stop = p.stop)
  return(df)
}

generate_oc_data <- function(posterior_probability_data, 
                             num_patients, 
                             confidence_go, 
                             confidence_stop) {
  
  # Find minimum responders that satisfy "Go!" criteria.
  responders_go <- posterior_probability_data[(
    posterior_probability_data$posterior_prob) > confidence_go,][1, 1]
  
  # Find minimum responders that satisfy "Stop!" criteria.
  responders_stop <- tail(posterior_probability_data[(
    posterior_probability_data$posterior_prob) <= confidence_stop,], n = 1)[1, 1]
  responders_stop <- ifelse(is.na(responders_stop), 0, responders_stop)
  
  # Operational characteristics (OC) probabilities.
  oc_probabilities <- lapply(0:100, function(x) oc_prob(x,
                                                        N = num_patients,
                                                        resp.go = responders_go,
                                                        resp.fut = responders_stop,
                                                        prob = 1-x/100))
  
  oc_data <- data.frame(do.call(rbind, oc_probabilities))
  
  return(oc_data)
}

create_oc_plot_data <- function(n.pat, 
                                MinPP,
                                conf.go,
                                conf.stop) {
  
  posterior_probability_data <- generate_posterior_probability_data(n.pat, MinPP) 
  
  oc_data <- generate_oc_data(posterior_probability_data, n.pat, conf.go, conf.stop) 
  
  # GO
  OCplot1 <- data.frame(
    "prob" = oc_data$prob,
    "prob.p" = oc_data$p.go,
    "col" = "green")
  OCplot1$group <- 1
  
  # EVAL
  OCplot2 <- data.frame(
    "prob"= oc_data$prob,
    "prob.p" = 1 - oc_data$p.stop - oc_data$p.go,
    "col" = "yellow")
  OCplot2$group <- 2
  
  # STOP
  OCplot3 <- data.frame(
    "prob" = oc_data$prob,
    "prob.p" = oc_data$p.stop,
    "col" = "red")
  OCplot3$group <- 3
  
  oc_plot_data <- rbind(OCplot1, OCplot2, OCplot3)
  oc_plot_data$group <- as.factor(oc_plot_data$group)
  
  return(oc_plot_data)
}

# Cumulative Decision Probabilities Plot.
draw_cumulative_decision_probabilities <- function(OCplota) {
  data <- OCplota[OCplota$prob<=1,]
  p <- ggplot(data=data) +
    geom_area(aes(x = prob, y = prob.p, fill = group), alpha=0.7) +
    scale_fill_manual(values = c('darkgreen', 'orange1', 'red2'),
                      labels = c("Go", "Evaluate", "Stop")) +  
    xlab("True Response Rate") + 
    ylab("Stop Probability") + 
    labs(fill = "Probabilities:") +  
    scale_x_continuous(expand=c(0,0),
                       breaks=seq(0,1,0.1),
                       labels=seq(0,1,0.1),
                       limits=c(-0.02,1.02)) + 
    scale_y_continuous(expand=c(0,0),
                       breaks=seq(0,1,0.1),
                       labels=seq(0,1,0.1),
                       limits=c(-0.09,1.09),
                       sec.axis = sec_axis(~.*1,breaks=seq(0,1,0.1),
                                           labels=seq(1,0,-.1),
                                           name="Go Probability")) + 
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 16),     
          legend.text = element_text(size = 14),  
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 17),
          plot.margin = margin(15, 15, 15, 15)) 
  return(p)
}

# Decision Probabilities Plot.
draw_decision_probabilities <- function(OCplota) {
  data <- OCplota[OCplota$prob<=1, ]
  p <- ggplot(data = data) + 
    geom_line(linewidth = 0.8,
              mapping = aes(x = prob, y = prob.p, col = group), alpha=0.4) +
    scale_color_manual("Probabilities:",
                       values = c("darkgreen",'orange1', 'red2'),
                       labels = c("Go","Evaluate","Stop")) +
    xlab("True Response Rate") + 
    ylab("Probability") + 
    scale_x_continuous(expand=c(0,0), 
                       breaks=seq(0,1,0.1), 
                       labels=seq(0,1,0.1),
                       limits=c(-0.02,1.02)) + 
    scale_y_continuous(expand=c(0,0), 
                       breaks=seq(0,1,0.1), 
                       labels=seq(0,1,0.1),
                       limits=c(-0.09,1.09)) + 
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 16),     
          legend.text = element_text(size = 14),  
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 17),
          plot.margin = margin(15, 15, 15, 15)) 
  return(p)
}

adapt_columns <- function(final_gating_table,
                          go_option,
                          eval_option,
                          stop_option) {
  colnames(final_gating_table) <- make.unique(colnames(final_gating_table))
  if (stop_option && go_option && !eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("EVALUATE", col)))))
  }
  if (!stop_option && go_option && eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("STOP", col))))) 
  }
  if (stop_option && !go_option && eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("GO", col)))))  
  }
  if (!stop_option && !go_option && eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("STOP", col))))) %>%
      dplyr::select(-which(sapply(., function(col) any(grepl("GO", col)))))
  }
  if (stop_option && !go_option && !eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("EVALUATE", col))))) %>%
      dplyr::select(-which(sapply(., function(col) any(grepl("GO", col)))))
  }
  if (!stop_option && go_option && !eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("STOP", col))))) %>%
      dplyr::select(-which(sapply(., function(col) any(grepl("EVALUATE", col)))))
  }
  if (!stop_option && !go_option && !eval_option) {
    final_gating_table <- final_gating_table %>% 
      dplyr::select(-which(sapply(., function(col) any(grepl("STOP", col))))) %>%
      dplyr::select(-which(sapply(., function(col) any(grepl("GO", col))))) %>%
      dplyr::select(-which(sapply(., function(col) any(grepl("EVALUATE", col)))))
  }
  return(final_gating_table)
}

make_decision_table_for_prob_of_interest <- function(n, 
                                                     tpp, 
                                                     confidence_go, 
                                                     confidence_stop, 
                                                     prob_of_interest,
                                                     stop_option,
                                                     go_option,
                                                     eval_option) {
  Nresp <- 0:n 
  orr_tpp_go <- tpp
  orr_tpp_stop <- tpp
  
  prob_min <- as.numeric(prob_of_interest[1])
  prob_max <- as.numeric(prob_of_interest[2])
  prob_of_interest <- seq(prob_min,prob_max,0.05)
  if (!(prob_max %in% prob_of_interest)) {
    prob_of_interest <- append(prob_of_interest, prob_max)
  }
  
  pp_go <- postprob(Nresp, # number of successes
                    n, # number of patients
                    orr_tpp_go, # threshold
                    parE = c(1, 1), # beta parameters, here uniform prior
                    log.p = FALSE) # log of the probability? (default: FALSE)
  
  min_responders_go <- min(Nresp[pp_go >= confidence_go])
  
  pp_stop <- postprob(Nresp, 
                      n, 
                      orr_tpp_stop, 
                      parE = c(1, 1), 
                      log.p = FALSE)
  
  max_responders_stop <- max(Nresp[pp_stop < confidence_stop])
  
  # distribution function for the binomial distribution with parameters size and prob.
  # stats::pbinom() is used to calculate the cumulative distribution function (CDF) 
  # values for a binomial distribution.
  binomial_prob_go <- stats::pbinom(min_responders_go - 1, # vector of quantiles
                                    n, # number of trials (zero or more).
                                    prob_of_interest, # probability of success on each trial
                                    lower.tail = FALSE, # TRUE: probabilities are P[X ≤ x], FALSE: P[X > x]
                                    log.p = FALSE)
  
  # distribution function for the binomial distribution with parameters size and prob.
  binomial_prob_stop <- stats::pbinom(max_responders_stop, 
                                      n, 
                                      prob_of_interest, 
                                      lower.tail = TRUE,  
                                      log.p = FALSE)
  
  binomial_prob_eval <- 1 - binomial_prob_go - binomial_prob_stop
  
  obs_orr_stop <- round(max_responders_stop/n*100,1)
  obs_orr_eval <- round(max_responders_stop/n*100,1)
  obs_orr_go <- round(min_responders_go/n*100,1)
  
  # Get number of responders for evaluate.
  eval_seq <- seq(max_responders_stop+1,min_responders_go-1,1)
  min_eval <- min(eval_seq)
  max_eval <- max(eval_seq)
  min_eval_obs_orr <- round(min_eval/n*100,1)
  max_eval_obs_orr <- round(max_eval/n*100,1)
  
  part1_first_col <- c('Decision Time Point',
                       'Number of Responders',
                       'Observed ORR',
                       'True ORR')
  part1_second_col <- c(paste0('N=',n),
                        paste0('0-',max_responders_stop),
                        paste0('0%-',obs_orr_stop,'%'),
                        'STOP')
  part1_third_col <- c(paste0('N=',n),
                       paste0(min_eval,'-',max_eval),
                       paste0(min_eval_obs_orr,'%-',max_eval_obs_orr,'%'),
                       'EVALUATE')
  part1_fourth_col <- c(paste0('N=',n),
                        paste0(min_responders_go,'-',n),
                        paste0(obs_orr_go,'%-100%'),
                        'GO')
  
  df1 <- data.frame(fcol=part1_first_col,scol=part1_second_col,tcol=part1_third_col,fcol=part1_fourth_col)
  df2 <- data.frame(fcol=paste0(percent_conditional_round_vectorized(prob_of_interest),'%'),
                    scol=paste0(percent_conditional_round_vectorized(binomial_prob_stop),'%'),
                    tcol=paste0(percent_conditional_round_vectorized(binomial_prob_eval),'%'),
                    fcol=paste0(percent_conditional_round_vectorized(binomial_prob_go),'%'))
  final_gating_table <- rbind(df1,df2)
  final_gating_table <- adapt_columns(final_gating_table,
                                      go_option,
                                      eval_option,
                                      stop_option)
  return(final_gating_table)
}

make_decision_table_for_multiple_sample_sizes <- function(n_vec, 
                                                          tpp, 
                                                          confidence_go, 
                                                          confidence_stop, 
                                                          prob_of_interest,
                                                          group_option,
                                                          stop_option,
                                                          go_option,
                                                          eval_option) {
  
  prob_min <- as.numeric(prob_of_interest[1])
  prob_max <- as.numeric(prob_of_interest[2])
  prob_of_interest <- seq(prob_min,prob_max,0.05)
  if (!(prob_max %in% prob_of_interest)) {
    prob_of_interest <- append(prob_of_interest, prob_max)
  }
  
  n_min <- n_vec[1]
  n_max <- n_vec[2]
  sample_sizes <- seq(n_min,n_max,5)
  if (!(n_max %in% sample_sizes)) {
    sample_sizes <- append(sample_sizes, n_max)
  }
  
  df_list <- list()
  go_df_list <- list()
  stop_df_list <- list()
  eval_df_list <- list()
  index <- 1
  for (n in sample_sizes) {
    Nresp <- 0:n 
    orr_tpp_go <- tpp
    orr_tpp_stop <- tpp
    pp_go <- postprob(Nresp, # number of successes
                      n, # number of patients
                      orr_tpp_go, # threshold
                      parE = c(1, 1), # beta parameters, here uniform prior
                      log.p = FALSE) # log of the probability? (default: FALSE)
    min_responders_go <- min(Nresp[pp_go >= confidence_go])
    pp_stop <- postprob(Nresp, 
                        n, 
                        orr_tpp_stop, 
                        parE = c(1, 1), 
                        log.p = FALSE)
    max_responders_stop <- max(Nresp[pp_stop < confidence_stop])
    binomial_prob_go <- stats::pbinom(min_responders_go - 1, # vector of quantiles
                                      n, # number of trials (zero or more).
                                      prob_of_interest, # probability of success on each trial
                                      lower.tail = FALSE, # TRUE: probabilities are P[X ≤ x], FALSE: P[X > x]
                                      log.p = FALSE)
    binomial_prob_stop <- stats::pbinom(max_responders_stop, 
                                        n, 
                                        prob_of_interest, 
                                        lower.tail = TRUE,  
                                        log.p = FALSE)
    binomial_prob_eval <- 1 - binomial_prob_go - binomial_prob_stop
    obs_orr_stop <- round(max_responders_stop/n*100,1)
    obs_orr_eval <- round(max_responders_stop/n*100,1)
    obs_orr_go <- round(min_responders_go/n*100,1)
    
    eval_seq <- seq(max_responders_stop+1,min_responders_go-1,1)
    min_eval <- min(eval_seq)
    max_eval <- max(eval_seq)
    min_eval_obs_orr <- round(min_eval/n*100,1)
    max_eval_obs_orr <- round(max_eval/n*100,1)
    
    part1_stop_col <- c(paste0('N=',n),
                        paste0('0-',max_responders_stop),
                        paste0('0%-',obs_orr_stop,'%'),
                        'STOP')
    
    part1_eval_col <- c(paste0('N=',n),
                        paste0(min_eval,'-',max_eval),
                        paste0(min_eval_obs_orr,'%-',max_eval_obs_orr,'%'),
                        'EVALUATE')
    
    part1_go_col <- c(paste0('N=',n),
                      paste0(min_responders_go,'-',n),
                      paste0(obs_orr_go,'%-100%'),
                      'GO')
    
    if (group_option) {
      stop_col <- c(part1_stop_col,paste0(percent_conditional_round_vectorized(binomial_prob_stop),'%'))
      eval_col <- c(part1_eval_col,paste0(percent_conditional_round_vectorized(binomial_prob_eval),'%'))
      go_col   <- c(part1_go_col,paste0(percent_conditional_round_vectorized(binomial_prob_go),'%'))
      df_stop <- data.frame(scol=stop_col)
      df_eval <- data.frame(tcol=eval_col)
      df_go <- data.frame(fcol=go_col)
      go_df_list[[index]] <- df_go
      stop_df_list[[index]] <- df_stop
      eval_df_list[[index]] <- df_eval
    } else {
      df1 <- data.frame(scol=part1_stop_col,tcol=part1_eval_col,fcol=part1_go_col)
      df2 <- data.frame(scol=paste0(percent_conditional_round_vectorized(binomial_prob_stop),'%'),
                        tcol=paste0(percent_conditional_round_vectorized(binomial_prob_eval),'%'),
                        fcol=paste0(percent_conditional_round_vectorized(binomial_prob_go),'%'))
      gating_table <- rbind(df1,df2)
      df_list[[index]] <- gating_table
    }  
    index <- index + 1
  }
  if (group_option) {
    combined_list <- lapply(list(stop_df_list, eval_df_list, go_df_list), function(df_list) {
      do.call(cbind, df_list)
    })
    combined_df <- do.call(cbind, combined_list)
    part1 <- c('Decision Time Point','Number of Responders','Observed ORR','True ORR')
    first_col <- c(part1, paste0(percent_conditional_round_vectorized(prob_of_interest),'%'))
    final_gating_table <- cbind(first_col, combined_df)
  } else {
    combined_df <- do.call(cbind, df_list)
    part1 <- c('Decision Time Point','Number of Responders','Observed ORR','True ORR')
    first_col <- c(part1, paste0(percent_conditional_round_vectorized(prob_of_interest),'%'))
    final_gating_table <- cbind(first_col, combined_df)
  }
  final_gating_table <- adapt_columns(final_gating_table,
                                      go_option,
                                      eval_option,
                                      stop_option)
  return(final_gating_table)
}

make_decision_table <- function(n, tpp, confidence_go, confidence_stop, prob_of_interest) {
  
  Nresp <- 0:n # Number of responders
  
  orr_tpp_go <- tpp
  orr_tpp_stop <- tpp
  
  min <- prob_of_interest[1]
  max <- prob_of_interest[2]
  prob_of_interest <- seq(min,max,0.05)
  
  pp_go <- postprob(Nresp, # number of successes
                    n, # number of patients
                    orr_tpp_go, # threshold
                    parE = c(1, 1), # beta parameters, here uniform prior
                    log.p = FALSE) # log of the probability? (default: FALSE)
  
  min_responders_go <- min(Nresp[pp_go >= confidence_go])
  
  pp_stop <- postprob(Nresp, 
                      n, 
                      orr_tpp_stop, 
                      parE = c(1, 1), 
                      log.p = FALSE)
  
  max_responders_stop <- max(Nresp[pp_stop < confidence_stop])
  
  # distribution function for the binomial distribution with parameters size and prob.
  # stats::pbinom() is used to calculate the cumulative distribution function (CDF) 
  # values for a binomial distribution.
  binomial_prob_go   <- stats::pbinom(min_responders_go - 1, # vector of quantiles
                                      n, # number of trials (zero or more).
                                      prob_of_interest, # probability of success on each trial
                                      lower.tail = FALSE, # TRUE: probabilities are P[X ≤ x], FALSE: P[X > x]
                                      log.p = FALSE)
  
  # distribution function for the binomial distribution with parameters size and prob.
  binomial_prob_stop <- stats::pbinom(max_responders_stop, 
                                      n, 
                                      prob_of_interest, 
                                      lower.tail = TRUE,  
                                      log.p = FALSE)
  
  binomial_prob_eval <- 1 - binomial_prob_go - binomial_prob_stop
  
  table_data = data.frame(
    N    = rep(n,length(prob_of_interest)),
    p    = prob_of_interest, 
    STOP = paste0(percent_conditional_round_vectorized(binomial_prob_stop),'%'), 
    EVAL = paste0(percent_conditional_round_vectorized(binomial_prob_eval),'%'), 
    GO   = paste0(percent_conditional_round_vectorized(binomial_prob_go),'%')
  )
  
  return(table_data)
}  
