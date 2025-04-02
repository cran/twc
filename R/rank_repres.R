#' Ensemble representativeness ranking
#'
#' The function \code{rank_repres} ranks the elements of a dataset ensemble in terms of a representativeness metric.
#'
#' @details
#' data is a data.table (time, value) method relates to the metric used to estimate the similarity to the ensemble mean/median. Can be mean, variance, correlation, slope, Kling–Gupta efficiency, Taylor skill score, Kullback–Leibler divergence, and all (for returning all of them)
#'
#' @import data.table
#' @importFrom stats approxfun cor density sd
#' @param data data.table with three variables: dataset, date, value
#' @param method character with seven options: mean, var (variance), slope, kge (Kling–Gupta Efficiency), tss (Taylor Skill Score), kld (Kullback–Leibler Divergence), all (default)
#' @param ensemble character with two options: mean, median (default)
#' @return data.table 2-column (dataset, repres_metric) or 8-column (dataset, mean, variance, ...) in case of method = all
#' @export
#' @examples
#' \dontrun{
#' download_data("gldas-vic", tempdir(), timestep = "yearly")
#' download_data("gpcc", tempdir(), timestep = "yearly")
#' download_data("ghcn", tempdir(), timestep = "yearly")
#' r1 <- raster::brick(paste0(tempdir(),
#' "/gldas-vic-v2-0_tp_mm_land_194801_201412_025_yearly.nc"))
#' s1 <- fldmean(r1)
#' s1$dataset <- "gldas-vic"
#' r2 <- raster::brick(paste0(tempdir(),
#' "/gpcc-v2022_tp_mm_land_198201_202012_025_yearly.nc"))
#' s2 <- fldmean(r2)
#' s2$dataset <- "gpcc"
#' r3 <- raster::brick(paste0(tempdir(),
#' "/ghcn-v2_tp_mm_land_190001_201505_025_yearly.nc"))
#' s3 <- fldmean(r2)
#' s3$dataset <- "ghcn"
#' dt <- rbind(r1, r2, r3)
#' dr <- rank_repres(dt[year(date) >= 1991 & year(date) <= 2012, .(dataset, date, value)])
#' }

rank_repres <- function(data, method = "all", ensemble = "median") {
  prec_ensemble <- data[, .(ensemble = match.fun(ensemble)(value, na.rm = TRUE)),
                        .(date)]

  if (method == "mean") {
    stat_ensemble <- prec_ensemble[, .(mean_ensemble = mean(ensemble, na.rm = TRUE))]
    prec_data <- data[, .(repres_metric = mean(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - mean_ensemble)/mean_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "var") {
    stat_ensemble <- prec_ensemble[, .(var_ensemble = sd(ensemble, na.rm = TRUE)^2)]
    prec_data <- data[, .(repres_metric = sd(value, na.rm = TRUE)^2), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - var_ensemble)/var_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "slope") {
    stat_ensemble <- prec_ensemble[, .(trend_ensemble = lm(ensemble ~ date)$coefficients[2])]
    prec_data <- data[, .(repres_metric = lm(value ~ date)$coefficients[2]), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - trend_ensemble)/trend_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "kge") {
    stat_ensemble <- prec_ensemble[, .(mean_ensemble = mean(ensemble, na.rm = TRUE),
                                       sd_ensemble = sd(ensemble, na.rm = TRUE))]
    stat_data <- data[, .(mean_prec = mean(value, na.rm = TRUE),
                          sd_prec = sd(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(stat_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(alpha = (sd_prec/mean_prec)/(sd_ensemble/mean_ensemble),
                                       beta = mean_prec/mean_ensemble), .(dataset)]
    prec_data <- merge(data, prec_ensemble, by = "date", allow.cartesian = TRUE)
    prec_data <- prec_data[, .(r_prec = cor(value, ensemble,
                                            use = "pairwise.complete.obs")),
                           .(dataset)]
    stat_ensemble <- merge(stat_ensemble, prec_data, by = "dataset")
    stat_ensemble <- stat_ensemble[, .(repres_metric = 1 - sqrt(((r_prec - 1)^2) + ((alpha - 1)^2) + ((beta - 1)^2))), .(dataset)]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "tss") {
    stat_ensemble <- prec_ensemble[, .(sd_ensemble = sd(ensemble, na.rm = TRUE))]
    stat_data <- data[, .(sd_prec = sd(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(stat_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(a = sd_prec/sd_ensemble), .(dataset)]
    prec_data <- merge(data, prec_ensemble, by = "date", allow.cartesian = TRUE)
    prec_data <- prec_data[, .(r_prec = cor(value, ensemble,
                                            use = "pairwise.complete.obs")),
                           .(dataset)]
    stat_ensemble <- merge(stat_ensemble, prec_data, by = "dataset")
    stat_ensemble <- stat_ensemble[, .(repres_metric = 2*(1 + r_prec)/((a + (1/a))^2)), .(dataset)]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "kld") {
    Q <- approxfun(density(prec_ensemble$ensemble))
    ALL_DATA <- unique(data$dataset)
    stat_ensemble <- foreach(idx = 1:length(ALL_DATA), .combine = rbind) %do% {
      prec_data <- data[dataset == ALL_DATA[idx]]
      P <- approxfun(density(prec_data$value))
      dx <- diff(sort(unique(prec_data$value)))
      dy <- diff(sort(unique(prec_ensemble$ensemble)))
      ex <- min(dx)
      ey <- min(dy)
      e <- min(ex, ey)/2
      n <- length(prec_data$value)
      x <- sort(prec_data$value)
      KL <- sum(log((P(x) - P(x - e))/(Q(x) - Q(x - e))), na.rm = TRUE) / n
      if (KL < 0) (KL <- 0)
      dummie <- data.table("dataset" = ALL_DATA[idx], repres_metric = KL)
      dummie
    }
    stat_ensemble[, repres_metric := 1 - repres_metric]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
  } else if (method == "all") {
    ## Mean
    stat_ensemble <- prec_ensemble[, .(mean_ensemble = mean(ensemble, na.rm = TRUE))]
    prec_data <- data[, .(repres_metric = mean(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - mean_ensemble)/mean_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    dummie_all <- stat_ensemble[order(-repres_metric)]
    dummie_all <- dummie_all[, .(mean = repres_metric), .(dataset)]
    ## Variance
    stat_ensemble <- prec_ensemble[, .(var_ensemble = sd(ensemble, na.rm = TRUE)^2)]
    prec_data <- data[, .(repres_metric = sd(value, na.rm = TRUE)^2), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - var_ensemble)/var_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
    stat_ensemble <- stat_ensemble[, .(dataset, variance = repres_metric)]
    dummie_all <- merge(dummie_all, stat_ensemble, by = "dataset")
    ## Slope
    stat_ensemble <- prec_ensemble[, .(trend_ensemble = lm(ensemble ~ date)$coefficients[2])]
    prec_data <- data[, .(repres_metric = lm(value ~ date)$coefficients[2]), .(dataset)]
    stat_ensemble <- cbind(prec_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(repres_metric = abs((repres_metric - trend_ensemble)/trend_ensemble)),
                                   .(dataset)]
    stat_ensemble <- stat_ensemble[, .(repres_metric = (1 - repres_metric)), .(dataset)]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
    stat_ensemble <- stat_ensemble[, .(dataset, slope = repres_metric)]
    dummie_all <- merge(dummie_all, stat_ensemble, by = "dataset")
    ## KGE
    stat_ensemble <- prec_ensemble[, .(mean_ensemble = mean(ensemble, na.rm = TRUE),
                                       sd_ensemble = sd(ensemble, na.rm = TRUE))]
    stat_data <- data[, .(mean_prec = mean(value, na.rm = TRUE),
                          sd_prec = sd(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(stat_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(alpha = (sd_prec/mean_prec)/(sd_ensemble/mean_ensemble),
                                       beta = mean_prec/mean_ensemble), .(dataset)]
    prec_data <- merge(data, prec_ensemble, by = "date", allow.cartesian = TRUE)
    prec_data <- prec_data[, .(r_prec = cor(value, ensemble,
                                            use = "pairwise.complete.obs")),
                           .(dataset)]
    stat_ensemble <- merge(stat_ensemble, prec_data, by = "dataset")
    stat_ensemble <- stat_ensemble[, .(repres_metric = 1 - sqrt(((r_prec - 1)^2) + ((alpha - 1)^2) + ((beta - 1)^2))), .(dataset)]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
    stat_ensemble <- stat_ensemble[, .(dataset, kge = repres_metric)]
    dummie_all <- merge(dummie_all, stat_ensemble, by = "dataset")
    ## TSS
    stat_ensemble <- prec_ensemble[, .(sd_ensemble = sd(ensemble, na.rm = TRUE))]
    stat_data <- data[, .(sd_prec = sd(value, na.rm = TRUE)), .(dataset)]
    stat_ensemble <- cbind(stat_data, stat_ensemble)
    stat_ensemble <- stat_ensemble[, .(a = sd_prec/sd_ensemble), .(dataset)]
    prec_data <- merge(data, prec_ensemble, by = "date", allow.cartesian = TRUE)
    prec_data <- prec_data[, .(r_prec = cor(value, ensemble,
                                            use = "pairwise.complete.obs")),
                           .(dataset)]
    stat_ensemble <- merge(stat_ensemble, prec_data, by = "dataset")
    stat_ensemble <- stat_ensemble[, .(repres_metric = 2*(1 + r_prec)/((a + (1/a))^2)), .(dataset)]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
    stat_ensemble <- stat_ensemble[, .(dataset, tss = repres_metric)]
    dummie_all <- merge(dummie_all, stat_ensemble, by = "dataset")
    ## KLD
    Q <- approxfun(density(prec_ensemble$ensemble))
    ALL_DATA <- unique(data$dataset)
    stat_ensemble <- foreach(idx = 1:length(ALL_DATA), .combine = rbind) %do% {
      prec_data <- data[dataset == ALL_DATA[idx]]
      P <- approxfun(density(prec_data$value))
      dx <- diff(sort(unique(prec_data$value)))
      dy <- diff(sort(unique(prec_ensemble$ensemble)))
      ex <- min(dx)
      ey <- min(dy)
      e <- min(ex, ey)/2
      n <- length(prec_data$value)
      x <- sort(prec_data$value)
      KL <- sum(log((P(x) - P(x - e))/(Q(x) - Q(x - e))), na.rm = TRUE) / n
      if (KL < 0) (KL <- 0)
      dummie <- data.table("dataset" = ALL_DATA[idx], repres_metric = KL)
      dummie
    }
    stat_ensemble[, repres_metric := 1 - repres_metric]
    stat_ensemble[repres_metric < 0, repres_metric := 0]
    stat_ensemble <- stat_ensemble[order(-repres_metric)]
    stat_ensemble <- stat_ensemble[, .(dataset, kld = repres_metric)]
    dummie_all <- merge(dummie_all, stat_ensemble, by = "dataset")
    stat_ensemble <- dummie_all
  }
  return(stat_ensemble)
}
