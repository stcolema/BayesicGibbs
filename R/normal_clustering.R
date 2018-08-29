#!/usr/bin/env Rscript

# === Functions ================================================================

# --- MCMC analysis ------------------------------------------------------------
#' @title entropy_window
#' @description  Find the point at which entropy stabilises in the iterations.
#'
#' @param entropy_vec A vector of numbers corresponding to entropy of each
#' iteration.
#' @param start An integer instructing which iteration to start with (default is
#' 1).
#' @param window_length The number of iterations to consider when considering
#' convergence (default is 25).
#' @param mean_tolerance A number. The threshold for how close the mean of the
#' two windows must be to be considered converged (default is 0.001).
#' @param sd_tolerance: A number. The threshold for how close the standard
#' deviation of the two windows must be to be considered converged (default is
#' 0.001).
#' @return The iteration at which convergence occurs in the clustering
entropy_window <- function(entropy_vec,
                           start = 1,
                           window_length = 25,
                           mean_tolerance = 0.001,
                           sd_tolerance = 0.001) {
  n <- length(entropy_vec)

  search_range <- seq(
    from = start,
    to = n - window_length,
    by = window_length
  )

  for (i in search_range) {

    # Create two windows looking forward from the current iteration and compare
    # their means and standard deviations
    win_1 <- entropy_vec[i:(i + window_length - 1)]
    win_2 <- entropy_vec[(i + window_length)
    :min((i + 2 * window_length - 1), n)]

    mean_1 <- mean(win_1)
    mean_2 <- mean(win_2)

    sd_1 <- sd(win_1)
    sd_2 <- sd(win_2)

    # If the differences are less than the predefined tolerances, return this
    # iteration as the point to burn up to
    if ((abs(mean_1 - mean_2) < mean_tolerance)
    + (abs(sd_1 - sd_2) < sd_tolerance)
    ) {
      return(i)
    }
  }
}

# --- Gibbs sampling -----------------------------------------------------------
#' Generates priors for the mean, degrees of freedom and scale parameters if not
#' set.
#'
#' @param data A matrix of the data being analysed.
#' @param mu_0 A d-vector. If NULL defaults to a vector of column means of data.
#' @param df_0 An integer. If NULL defaults to d + 2.
#' @param scale_0 A positive definite matrix. The prior for the scale parameter
#' of the inverse wishart distribution; if NULL defaults to a diagonal matrix.
#' @param N The number of entries in data.
#' @param k The number of clusters used.
#' @param d The number of columns in data.
#' @return A named list of the three hyperparameters, mean, scale and degrees of
#'  freedom
empirical_bayes_initialise <- function(data, mu_0, df_0, scale_0, N, k, d) {
  parameters <- list()
  if (is.null(mu_0)) {
    mu_0 <- colMeans(data)
  }

  if (is.null(df_0)) {
    df_0 <- d + 2
  }

  if (is.null(scale_0)) {
    scale_0 <- diag(colSums((data - mean(data))^2) / N) / (k^(1 / d))
    if (any(is.na(scale_0))) {
      scale_0 <- diag(d) / (k^(1 / d))
    }
  }
  parameters$mu_0 <- mu_0
  parameters$df_0 <- df_0
  parameters$scale_0 <- scale_0

  return(parameters)
}

#' @title Gibbs sampling
#' @description Carries out gibbs sampling of data and returns a similarity matrix for points
#'
#' @param data A matrix of the data being analysed.
#' @param k The number of clusters.
#' @param class_labels A vector of unsigned integers representing the initial
#' cluster of the corresponding point in data
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param mu_0 A d-vector; prior on mean. If NULL defaults to mean of data.
#' @param df_0 The prior on the degrees of freedom. if NULL defaults to d + 2.
#' @param scale_0 The prior on the scale for the Inverse Wishart. If NULL
#' generated using an empirical method.
#' @param lambda_0 The prior of shrinkage for mean distribution.
#' @param concentration_0 The prior for dirichlet distribution of cluster
#' weights.
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param outlier A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution
#' @param t_df The degrees of freedom for the outlier t-distribution (default
#' is 4)
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is FALSE)
gibbs_sampling <- function(data, k, class_labels, fix_vec,
                           d = NULL,
                           N = NULL,
                           num_iter = NULL,
                           burn = 0,
                           mu_0 = NULL,
                           df_0 = NULL,
                           scale_0 = NULL,
                           lambda_0 = 0.01,
                           concentration_0 = 0.1,
                           thinning = 25,
                           outlier = FALSE,
                           t_df = 4.0,
                           record_posteriors = FALSE) {
  if (is.null(d)) {
    d <- ncol(data)
  }

  if (is.null(N)) {
    N <- nrow(data)
  }


  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  if (thinning > (num_iter - burn)) {
    if (thinning > (num_iter - burn) & thinning < 5 * (num_iter - burn)) {
      stop("Thinning factor exceeds iterations feasibly recorded. Stopping.")
    } else if (thinning > 5 * (num_iter - burn) & thinning < 10 * (num_iter - burn)) {
      stop("Thinning factor relatively large to effective iterations. Stopping algorithm.")
    } else {
      warning(paste0(
        "Thinning factor relatively large to effective iterations.",
        "\nSome samples recorded. Continuing but please check input"
      ))
    }
  }

  data <- as.matrix(data)

  # Empirical Bayes
  parameters_0 <- empirical_bayes_initialise(data, mu_0, df_0, scale_0, N, k, d)

  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0

  if (is.null(concentration_0)) {
    concentration_0 <- rep(0.1, k)
  } else if (length(concentration_0) < k) {
    print(paste0(
      "Creating vector of ", k + outlier, " repetitions of ", concentration_0,
      "Creating vector of ", k + outlier, " repetitions of ", concentration_0,
      " for concentration prior."
    ))
    concentration_0 <- rep(concentration_0, k + outlier)
  }

  sim <- gaussian_clustering(
    num_iter,
    concentration_0,
    scale_0,
    class_labels,
    fix_vec,
    mu_0,
    lambda_0,
    data,
    df_0,
    k,
    burn,
    thinning,
    outlier,
    t_df,
    record_posteriors
  )
}

# --- Categorical clustering ---------------------------------------------------
#' @title Phi prior
#' @description Generates a prior for the phi vector for each variable for the Dirichlet
#' distribution
#' @param matrix_data A matrix of data.
#' @return A list of vectors of the proportion of each level across all of
#' matrix_data.
phi_prior <- function(matrix_data) {

  # lambda function applies ``table'' to each column of matr_data before
  # removing names. This ensures the same output type for the case when all
  # variables have the same number of levels and when they do not (rather than
  # merely using one or two lapply's of table)
  unnamed_list_prop <- lapply(
    1:ncol(matrix_data),
    function(i) {
      table(matrix_data[, i])
    }
  ) %>%
    unname() %>%
    lapply("/", nrow(matrix_data))
}

# --- Heatmap ------------------------------------------------------------------
#' @title Annotated Heatmap
#' @description Returns and prints an annotated pheatmap
#'
#' @param input_data A matrix of data to be heat mapped. Needs column and
#' rownames.
#' @param annotation_row A data frame of the annotation variable(s). Names must
#' match with input_data. If NULL returns normal pheatmap.
#' @param sort_by_col The name of the column to sort input_data and
#' annotation_row by when heatmapping. If pheatmap is instructed to sort_rows
#' this has no impact. Default is NULL in which case no sorting occurs.
#' @param train If FALSE returns normal pheatmap.
#' @param ... The usual inputs for pheatmap.
#' @return An annotated pheatmap from the pheatmap package
annotated_heatmap <- function(input_data, annotation_row = NULL,
                              sort_by_col = NULL,
                              train = NULL,
                              ...) {
  if (is.null(annotation_row) & !(isTRUE(train) | is.null(train))) {
    stop("If data")
  }
  dissim <- input_data

  sort_by_col <- attempt::try_catch(
    expr = enquo(sort_by_col),
    .e = NULL,
    .w = NULL
  )

  # print(sort_by_col)

  if (sort_by_col != quo(NULL)) {
    # print("HIII")

    # sort_col <- enquo(sort_by_col)

    col_of_interest <- annotation_row %>%
      dplyr::select(!!sort_by_col)

    if (!is.null(annotation_row)) {
      combined_data <- bind_cols(dissim, annotation_row)

      sorted_data <- combined_data %>%
        arrange(!!sort_by_col)

      annotation_row <- sorted_data %>%
        dplyr::select(one_of(names(annotation_row)))
      rownames(annotation_row) <- rownames(input_data)

      dissim <- sorted_data %>%
        dplyr::select(-one_of(names(annotation_row)))
    } else {
      dissim <- dissim %>%
        arrange(!!sort_by_col)
    }

    rownames(dissim) <- rownames(input_data)
  }

  # Colour scheme for heatmap
  col_pal <- RColorBrewer::brewer.pal(9, "Blues")

  if (!is.null(annotation_row) | !(is.null(train) | isTRUE(train))) {
    feature_names <- names(annotation_row)

    # Annotation colours
    new_cols_list <- list()
    my_colours <- list()
    for (feature in feature_names) {
      outlier_present <- FALSE

      # print(feature)
      types_feature_present <- unique(annotation_row[[feature]][!is.na(annotation_row[[feature]])])

      # print(types_feature_present)

      if (feature == "Predicted_class") {
        if ("Outlier" %in% types_feature_present) {
          outlier_present <- TRUE
          types_feature_present <- types_feature_present[types_feature_present != "Outlier"]
        }
      }

      # features_present[[feature]] <- types_feature_present
      new_cols_list[[feature]] <- colorRampPalette(grDevices::rainbow(length(types_feature_present)))



      my_colours[[feature]] <- new_cols_list[[feature]](length(types_feature_present))

      if (outlier_present) {
        my_colours[[feature]] <- c(my_colours[[feature]], "black")
        types_feature_present <- c(types_feature_present, "Outlier")
      }

      names(my_colours[[feature]]) <- types_feature_present
    }

    # Heatmap
    if (is.null(train) | isTRUE(train)) {
      heat_map <- pheatmap(dissim,
        annotation_row = annotation_row,
        annotation_colors = my_colours,
        ...
      )
    }
  }
  else {
    heat_map <- pheatmap(
      dissim,
      ...
    )
  }
  return(heat_map)
}

# === Wrapper ==================================================================
#' @title MCMC out
#' @description Returns mean, variance and similarity posteriors from Gibbs sampling with
#' option of pheatmap
#'
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param class_labels_0 An optional prior for clusters in MS_object. If NULL
#' defaults to a randomly generated set using the proportion in the labelled
#' data.
#' @param mu_0 A d-vector; prior on mean. If NULL defaults to mean of data.
#' @param df_0 The prior on the degrees of freedom. if NULL defaults to d + 2.
#' @param scale_0 The prior on the scale for the Inverse Wishart. If NULL
#' generated using an empirical method.
#' @param lambda_0 The prior of shrinkage for mean distribution.
#' @param concentration_0 The prior for dirichlet distribution of cluster
#' weights.
#' @param train: instruction to include all data (NULL), labelled data (TRUE) or
#' unlabelled data (FALSE). Default is NULL.
#' @param num_iter The number of iterations to sample over.
#' @param burn The number of iterations to record after (i.e. the burn-in).
#' @param thinning The step between iterations for which results are recorded in
#' the mcmc output.
#' @param heat_plot A bool. Instructs saving and printing of heatmap of
#' similarity matrix from Gibbs sampling. Default is TRUE.
#' @param main String. The title for heatmap, default is "heatmap_for_similarity".
#' @param cluster_row A bool. Instructs pheatmap to cluster rows using a tree.
#' @param cluster_cols: A bool. instructs pheatmap to cluster columns using a
#' tree.
#' @param fontsize: The size of font in pheatmap.
#' @param fontsize_row: The fontsize in rows in pheatmap.
#' @param fontsize_col: The fontsize in columns in pheatmap.
#' @param entropy_plot A bool instructing function to save a plot of the entropy
#' across all iterations. Default is TRUE.
#' @param window_length A number. Input to entropy ploy function. Default is
#' min(25, num_iter / 5).
#' @param mean_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param sd_tolerance Input to entropy ploy function. Default is 0.0005.
#' @param outlier A bool instructing the sampler to consider an additional
#' outlier cluster following a t-distribution
#' @param t_df The degrees of freedom for the outlier t-distribution (default
#' is 4)
#' @param prediction_threshold The minimum proportion of recorded iterations
#' for which a point is in its most common cluster for which a prediction is
#' returned. If below this predicted class is NA.
#' @param record_posteriors A bool instructing the mcmc function to record the
#' posterior distributions of the mean and variance for each cluster
#' (default is FALSE)
#' @return A named list including at least the output from the gibbs sampler,
#' but can include two pheatmaps and a scatter plot of the entropy over
#' iterations.
#' @examples
#' data("hyperLOPIT2015") # MS object from pRolocData
#' mcmc_object <- mcmc_out(hyperLOPIT2015,
#'   num_iter = 10000,
#'   burn = 1000,
#'   thinning = 50,
#'   outlier = TRUE,
#'   heat_plot = TRUE,
#'   main = "Gene clustering by organelle",
#'   prediction_threshold = 0.5
#' )
mcmc_out <- function(MS_object,
                     class_labels_0 = NULL,
                     mu_0 = NULL,
                     df_0 = NULL,
                     scale_0 = NULL,
                     lambda_0 = 0.01,
                     concentration_0 = 0.1,
                     train = NULL,
                     num_iter = NULL,
                     burn = NULL,
                     thinning = 25,
                     heat_plot = TRUE,
                     main = "heatmap_for_similarity",
                     cluster_row = T,
                     cluster_cols = T,
                     fontsize = 10,
                     fontsize_row = 6,
                     fontsize_col = 6,
                     entropy_plot = TRUE,
                     window_length = min(25, num_iter / 5),
                     mean_tolerance = 0.0005,
                     sd_tolerance = 0.0005,
                     outlier = FALSE,
                     t_df = 4.0,
                     prediction_threshold = 0.6,
                     record_posteriors = FALSE) {
  # MS data
  MS_data <- MS_dataset(MS_object, train = train)

  mydata <- MS_data$data
  nk <- MS_data$nk
  row_names <- MS_data$row_names
  fix_vec <- MS_data$fix_vec

  class_labels <- data.frame(Class = mydata$markers)

  classes_present <- unique(fData(markerMSnSet(MS_object))[, "markers"])

  rownames(class_labels) <- rownames(mydata)

  # Numerical data of interest for clustering
  num_data <- mydata %>%
    dplyr::select(-markers)

  # Parameters
  k <- length(classes_present)
  N <- nrow(num_data)
  d <- ncol(num_data)

  # Key to transforming from int to class
  class_labels_key <- data.frame(Class = classes_present) # , Class_num = 1:k)
  class_labels_key %<>%
    arrange(Class) %>%
    dplyr::mutate(Class_key = as.numeric(Class))

  class_labels %<>%
    mutate(Class_ind = as.numeric(mydata$markers))

  if (outlier) {
    outlier_row <- data.frame(Class = c("Outlier"), Class_key = c(k + 1))
    class_labels_key <- suppressWarnings(bind_rows(class_labels_key, outlier_row))
  }

  class_labels_0 <- cluster_label_prior(class_labels_0, nk, train, MS_object, N)


  if (is.null(num_iter)) {
    num_iter <- floor(min((d^2) * 1000 / sqrt(N), 10000))
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  if (thinning > (num_iter - burn)) {
    if (thinning > (num_iter - burn) & thinning < 5 * (num_iter - burn)) {
      stop("Thinning factor exceeds iterations feasibly recorded. Stopping.")
    } else if (thinning > 5 * (num_iter - burn) & thinning < 10 * (num_iter - burn)) {
      stop("Thinning factor relatively large to effective iterations. Stopping algorithm.")
    } else {
      warning(paste0(
        "Thinning factor relatively large to effective iterations.",
        "\nSome samples recorded. Continuing but please check input"
      ))
    }
  }

  gibbs <- gibbs_sampling(num_data, k, class_labels_0, fix_vec,
    d = d,
    N = N,
    num_iter = num_iter,
    burn = burn,
    mu_0 = mu_0,
    df_0 = df_0,
    scale_0 = scale_0,
    lambda_0 = lambda_0,
    concentration_0 = concentration_0,
    thinning = thinning,
    outlier = outlier,
    t_df = t_df,
    record_posteriors = record_posteriors
  )

  print("Gibbs sampling complete")


  # Create a dataframe for the predicted class
  class_allocation_table <- with(
    stack(data.frame(t(gibbs$class_record))),
    table(ind, values)
  )

  eff_iter <- ceiling((num_iter - burn) / thinning)

  # Create a column Class_key containing an integer in 1:k representing the most
  # common class allocation, and a Count column with the proportion of times the
  # entry was allocated to said class
  predicted_classes <- data.frame(
    Class_key =
      as.numeric(colnames(class_allocation_table)
      [apply(
          class_allocation_table,
          1,
          which.max
        )]),
    Count = apply(class_allocation_table, 1, max) / eff_iter
  )

  # Change the prediction to NA for any entry with a proportion below the input
  # threshold
  predicted_classes[predicted_classes$Count < prediction_threshold, ] <- NA

  predicted_classes$Class <- class_labels_key$Class[match(
    predicted_classes$Class_key,
    class_labels_key$Class_key
  )]

  gibbs$predicted_class <- predicted_classes

  # Example input for annotation_row in pheatmap
  annotation_row <- class_labels %>% dplyr::select(Class)
  annotation_row %<>%
    mutate(Predicted_class = predicted_classes$Class)

  rownames(num_data) <- row_names
  # print(rownames(mydata[1:10,]))

  rownames(annotation_row) <- rownames(num_data)

  col_pal <- RColorBrewer::brewer.pal(9, "Blues")

  pauls_heatmap <- annotated_heatmap(num_data, annotation_row,
    sort_by_col = Predicted_class,
    train = train,
    main = "Paul's sense check heatmap",
    cluster_row = FALSE,
    cluster_cols = FALSE,
    color = col_pal,
    fontsize = fontsize,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col
  )

  # return(pauls_heatmap)

  all_data <- dplyr::bind_cols(num_data, dplyr::select(gibbs$predicted_class, Class))
  all_data$Fixed <- fix_vec
  all_data$Protein <- rownames(num_data)
  # rownames(all_data) <- rownames(num_data)

  if (heat_plot) {

    # dissimilarity matrix
    dissim <- 1 - gibbs$similarity

    # Require names to associate data in annotation columns with original data
    colnames(dissim) <- rownames(num_data)
    rownames(dissim) <- rownames(num_data)

    col_pal <- RColorBrewer::brewer.pal(9, "Blues")

    heat_map <- annotated_heatmap(dissim, annotation_row,
      train = train,
      main = main,
      cluster_row = cluster_row,
      cluster_cols = cluster_cols,
      color = col_pal,
      fontsize = fontsize,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col
    )
  }
  if (entropy_plot) {
    entropy_data <- data.frame(Index = 1:num_iter, Entropy = gibbs$entropy)

    rec_burn <- entropy_window(gibbs$entropy,
      window_length = window_length,
      mean_tolerance = mean_tolerance,
      sd_tolerance = sd_tolerance
    )

    # Check if instantly ok
    rec_burn <- ifelse(is.null(rec_burn), 1, rec_burn)

    entropy_scatter <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
      geom_point() +
      geom_vline(mapping = aes(xintercept = rec_burn, colour = "Reccomended"), lty = 2) +
      geom_vline(mapping = aes(xintercept = burn, colour = "Implemented"), lty = 4) +
      ggtitle("Entropy over iterations including recommended and implemented burn") +
      xlab("Iteration") + ylab("Entropy") +
      scale_color_manual(name = "Burn", values = c(
        Reccomended = "red",
        Implemented = "blue"
      ))
  }
  if (heat_plot & entropy_plot) {
    return(list(
      gibbs = gibbs,
      heat_map = heat_map,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
    ))
  }
  if (heat_plot) {
    return(list(
      gibbs = gibbs,
      heatmap = heat_map,
      data = all_data
    ))
  }
  if (entropy_plot) {
    return(list(
      gibbs = gibbs,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn,
      data = all_data
    ))
  }

  return(list(gibbs = gibbs, data = all_data))
}

#' @title MS dataset
#' @description Given an MS object from pRolocData, returns the numerical data,
#' the count of each class and the vector recording which labels are fixed (i.e
#' from the training set) and unfixed.
#'
#' @param MS_object A dataset in the format used by pRolocdata.
#' @param train: instruction to include all data (NULL), labelled data (TRUE) or
#' unlabelled data (FALSE). Default is NULL.
#' @examples
#' data("hyperLOPIT2015") # MS object from pRolocData
#' MS_data <- MS_dataset(hyperLOPIT2015)
MS_dataset <- function(MS_object, train = NULL) {
  # Data with labels
  mydata_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = TRUE
  )

  fixed <- rep(TRUE, nrow(mydata_labels))

  mydata_no_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = FALSE
  )

  not_fixed <- rep(FALSE, nrow(mydata_no_labels))

  nk <- tabulate(fData(markerMSnSet(MS_object))[, "markers"])

  mydata_no_labels$markers <- NA

  if (is.null(train)) {
    row_names <- c(rownames(mydata_labels), rownames(mydata_no_labels))
    mydata <- suppressWarnings(bind_rows(mydata_labels, mydata_no_labels))
    fix_vec <- c(fixed, not_fixed)
  } else if (isTRUE(train)) {
    row_names <- c(rownames(mydata_labels))
    mydata <- mydata_labels
    fix_vec <- fixed
  } else {
    row_names <- c(rownames(mydata_no_labels))
    mydata <- mydata_no_labels
    fix_vec <- not_fixed
  }
  return(list(data = mydata, fix_vec = fix_vec, row_names = row_names, nk = nk))
}



#' @title Cluster label prior
#' @description Generates a vector of labels if required
#' @param class_labels_0 An optional prior for clusters in MS_object. If NULL
#' defaults to a randomly generated set using the proportion in the labelled
#' data.
#' @param nk Table of the frequency of the classes in the datset
#' @param train A NULL, TRUE or FALSE value informing if supervised (TRUE),
#' semi-supervised (NULL) or unsupervised (FALSE)
#' @param MS_object The MS object from pRolocData being analysed
#' @param k The number of clusters in the dataset
#' @param N The number of observations
#' @return A vector of integers corresponding to the cluster allocation of the N
#' observations
cluster_label_prior <- function(class_labels_0,
                                nk,
                                train,
                                MS_object,
                                k,
                                N) {
  # Generate class labels
  if (is.null(class_labels_0)) {
    class_weights <- nk / sum(nk)
    if (is.null(train)) {
      fixed_labels <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
      class_labels_0 <- c(fixed_labels, sample(1:k, N - length(fixed_labels),
        replace = T,
        prob = class_weights
      ))


      # class_labels_0 <- c(fixed_labels, sample(1:k, nrow(mydata_no_labels),
      #                                          replace = T,
      #                                          prob = class_weights
      # ))
    } else if (isTRUE(train)) {
      class_labels_0 <- as.numeric(fData(markerMSnSet(MS_object))[, "markers"])
    } else {
      class_labels_0 <- sample(1:k, N, replace = T, prob = class_weights)
    }
  }
  return(class_labels_0)
}

# === Olly =====================================================================
#
# set.seed(5)
#
# # MS object
# data("HEK293T2011") # Human Embroyonic Kidney dataset
# data("hyperLOPIT2015") # Olly's normal data I think
# t1 <- Sys.time()
#
# stuff <- mcmc_out(HEK293T2011,
#   num_iter = 10000,
#   burn = 1000,
#   thinning = 50,
#   outlier = TRUE,
#   heat_plot = TRUE,
#   main = "Gene clustering by organelle",
#   prediction_threshold = 0.5
# )
#
# t2 <- Sys.time()
#
# t2 - t1 # how long does it take
#
# # To plot the entropy over iterations
# # stuff$entropy_plot
# # str(stuff$gibbs$class_prob)
# #
# # str(stuff$gibbs$class_record)
# #
# # stuff$gibbs$predicted_class
# #
# # y <- stuff$gibbs$class_record
# # z <- with(stack(data.frame(t(y))), table(ind, values))
# # predicted_classes <- data.frame(Class_key = as.numeric(colnames(z)[apply(z, 1, which.max)]),
# #                                 Count = apply(z, 1, max))
# #
# # predicted_classes$Class <- df2$B[match(df1$Class_key, df2$Class_key)]
# #
# # predicted_classes %<>%
# #   dplyr::mutate(Class = 0)
# #
# # summary(stuff$gibbs$predicted_class)
#
#
# # stuff$Protein <- rownames(stuff)
# z <- melt(stuff$data, id.vars = c("Protein", "Class", "Fixed"))
# # z$Fixed <- as.numeric(z$Fixed)
#
# curr_data <- z %>%
#   filter(Class %in% c("Cytosol", "PM"))
#
# curr_data %<>% group_by(Protein)
#
# ggplot(data = curr_data,
#        aes(x = variable, y = value, colour = Class, group = Protein)
#        ) +
#   geom_point() +
#   geom_line(aes(linetype = Fixed)) +
#   # gghighlight(! Fixed, use_direct_label = FALSE) +
#   NULL
#
