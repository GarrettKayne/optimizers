next_trading_day <- function(x) {

  tbl <- GCAMCQT::read_gildata("qt_tradingdaynew", .load = FALSE)[["qt_tradingdaynew"]]
  tbl <- tbl[J(90)][IFTRADINGDAY == 1, .(DATE = TRADINGDATE, NEW = TRADINGDATE)]
  data.table::setkey(tbl, DATE)
  tbl[J(x + 1), NEW, roll = -Inf]
}


size_factor <- function(mv) {

  ln_mv <- log(mv)
  mean_mv <- weighted.mean(ln_mv, mv)
  sd_mv <- sd(ln_mv)
  (ln_mv - mean_mv) / sd_mv
}


add_size_factor <- function(model) {

  model[, SIZE := size_factor(MV_Float), keyby = EndDate]
}



pct_rank <-  function(x) {

  stopifnot(is.numeric(x))
  min_pct <- 0.30
  total_size <- length(x)
  eff_size <- sum(!is.na(x))

  if (eff_size == 1 || eff_size / total_size < min_pct) {
    res <- rep(NA_real_, length(x))
  } else {
    res <- (frank(-x, na.last = "keep", ties.method = "average") - 1) / (eff_size - 1)
  }
  res
}


normalize_rk <- function(rk) {

  pct_rk <- pct_rank(rk)
  lower <- pnorm(-3)
  upper <- 1 - lower
  scaled_rk <- (upper - lower) * pct_rk + lower # scale it to 0.01 ~ 0.99
  qnorm(1 - scaled_rk)
}


normalize_model_score <- function(model) {

  model[, SCORE_RK := pct_rank(ModelScore), keyby = EndDate]
  model[, SCORE := normalize_rk(ModelScore), keyby = EndDate]
}


gildata_index_code <- function(bmk) {

  switch(
    bmk,
    CSI300 = 3145L,
    CSI500 = 4978L,
    CSI800 = 4982L,
    stop("Undefined bchmk.", .call = FALSE)
  )
}


bmk_size <- function(bchmk = c("CSI300", "CSI500", "CSI800"), model) {
  bchmk <- match.arg(bchmk)
  index <- switch(
    bchmk,
    CSI300 = 3145L,
    CSI500 = 4978L,
    CSI800 = 4982L,
    stop("Undefined bchmk.", .call = FALSE)
  )

  dates <- sort(unique(model$EndDate))
  indexUniv <- sparkling:::f_index_univ(index_inner_code = index, dates)
  readData(c("jy_icmpwt"), "dataPool")
  indexUniv <- jy_icmpwt[IndexCode == index, .(InnerCode, EndDate, Weight)][indexUniv, on = c("InnerCode", "EndDate"), roll = TRUE]
  indexUniv <- indexUniv[!is.na(Weight)]
  setnames(indexUniv, c("INNER_CODE", "DATE", "WEIGHT"))

  size <- model[indexUniv, SIZE, on = c("InnerCode" = "INNER_CODE", "EndDate" = "DATE"), roll = TRUE]
  indexUniv[, SIZE := size]
  indexUniv <- indexUniv[!is.na(SIZE)]
  indexUniv[, WEIGHT := WEIGHT / sum(WEIGHT), by = DATE]
  list(
    size = indexUniv[, .(SIZE = weighted.mean(SIZE, WEIGHT, na.rm = TRUE)), keyby = DATE],
    raw = indexUniv
  )
}

cdar_optimizer_frapo <- function(ini_weight, obj,
                                 group, cum_rtn,
                                 cvar_alpha, cvar_upper,
                                 group_bound,
                                 max_sec_weight,
                                 max_total_weight,
                               turnover_bound) {
  cum_rtn <- as.data.frame(cum_rtn)
  cum_rtn1 <- cum_rtn[2:nrow(cum_rtn),]
  cum_rtn[2:nrow(cum_rtn),] <- cum_rtn1-cum_rtn[1:nrow(cum_rtn)-1,]
  cum_rtn <- as.data.table(cum_rtn)

  no_obs <- nrow(cum_rtn)
  no_security <- length(obj)
  # stopifnot(no_security >= 1, no_obs >= 2, isTRUE(all.equal(names(obj), colnames(cum_rtn))))
  security_names <- names(obj)
  cum_rtn <- unname(zoo::coredata(cum_rtn))


  Opti_obj <-
    cbind(0,
          matrix(obj, 1, no_security),
          matrix(0, 1, no_obs),
          # matrix(0, 1, no_obs),
          matrix(0, 1, no_security))
  c1 <- cbind(
    1,
    matrix(0, 1, no_security),
    1 / ((1 - cvar_alpha) * no_obs) * matrix(1, 1, no_obs),
    # matrix(0, 1, no_obs),
    matrix(0, 1, no_security)
  )
  c2 <-
    cbind(-matrix(1, no_obs, 1),
          -as.matrix(cum_rtn),
          -diag(rep(1, no_obs)),
          # diag(rep(1, no_obs)),
          matrix(0, no_obs, no_security))
  c5 <- cbind(matrix(0,no_obs,1),
              matrix(0, no_obs, no_security),
              -diag(rep(1, no_obs)),
              # matrix(0, no_obs, no_obs),
              matrix(0, no_obs, no_security)
  )
  c6 <-
    cbind(0,
          matrix(1, 1, no_security),
          matrix(0, 1, no_obs),
          # matrix(0, 1, no_obs),
          matrix(0, 1, no_security))
  c7 <- cbind(0,
              matrix(0, 1, no_security),
              matrix(0, 1, no_obs),
              # matrix(0, 1, no_obs),
              matrix(1, 1, no_security))
  c8 <- cbind(matrix(0,no_security,1),
              diag(rep(1, no_security)),
              matrix(0, no_security, no_obs),
              # matrix(0, no_security, no_obs),
              -diag(rep(1, no_security))
  )
  c9 <- cbind(matrix(0,no_security,1),
              -diag(rep(1, no_security)),
              matrix(0, no_security, no_obs),
              # matrix(0, no_security, no_obs),
              -diag(rep(1, no_security))
  )
  c10 <- cbind(matrix(0,no_security,1),
               matrix(0, no_security, no_security),
               matrix(0, no_security, no_obs),
               # matrix(0, no_security, no_obs),
               -diag(rep(1, no_security))
  )

  # The constraint of Industry
  H <- matrix(0, length(group), 1 + 2*no_security + no_obs)
  for (i in 1:length(group)) {
    for (j in 1:length(security_names)) {
      if (as.numeric(security_names[j]) %in% group[[i]]) {
        H[i, j + 1] <- 1
      }
    }
  }

  mat <- rbind(c1, c2, c5, c6, H, c7, c8, c9, c10)
  dir <- matrix("<=", 3*no_security+2 * no_obs + 3 + length(group), 1)

  rhs <-
    cbind(
      cvar_upper,
      matrix(0, 1, 2 * no_obs),
      max_total_weight,
      matrix(group_bound, 1, length(group)),
      turnover_bound,
      matrix(ini_weight, 1, no_security),
      -matrix(ini_weight, 1, no_security),
      matrix(0, 1, no_security)
    )
  bounds <- list(lower = list(ind = c(1L), val = c(-Inf)),
                 upper = list(ind = c(2:(no_security + 1)), val = rep(max_sec_weight, no_security)))

  opt_res <-
    Rglpk::Rglpk_solve_LP(Opti_obj, mat, dir, rhs, bounds, max = TRUE)

  weights <- round(opt_res$solution[2:(no_security+1)], 6L)
  names(weights) <- security_names

  cvar <- local({
    # weights <- opt_res$solution[2:(no_security+1)]
    # equity <- matrix(apply(cum_rtn, 1, function(x) sum(x * weights)), ncol = 1)
    # uvals <- opt_res$solution[(no_security + no_obs + 2):(no_security + 2*no_obs +1)]
    # dd <- uvals - equity
    opt_res$solution[1]*c1
    # mean(dd[dd >= z])
  })

  list(
    weight =
      data.table(
        INNER_CODE = as.integer(names(weights)[weights > 0]),
        WEIGHT = weights[weights > 0.0],
        key = "INNER_CODE"
      ),
    cvar = cvar
  )

}


prev_trading_day <- function(x, n) {

  tbl <- GCAMCQT::read_gildata("qt_tradingdaynew", .load = FALSE)[["qt_tradingdaynew"]]
  tbl <- tbl[J(90)][IFTRADINGDAY == 1, .(DATE = TRADINGDATE, NEW = TRADINGDATE)]
  data.table::setkey(tbl, DATE)
  t_1 <- tbl[J(x), NEW, roll = Inf]
  tbl[which(tbl$DATE == t_1) - n, NEW]
}



cal_cum_rtn <- function(candle, date, inner_codes) {

  date_from <- prev_trading_day(date - 1, 220L)
  rtn <- candle[DATE >= date_from & DATE < date & INNER_CODE %in% inner_codes,
                .(DATE, INNER_CODE, CLOSE_ADJ = CLOSE * ADJ_FACTOR)]
  rtn[, c("CUM_RTN", "CT", "ZERO_CT") := .(cumr <- CLOSE_ADJ / CLOSE_ADJ[1L] - 1.0, .N, sum(round(cumr, 2) == 0.0)),
      by = INNER_CODE]
  max_zero_pct <- 0.2
  # browser()
  # cum_rtn <- as.data.frame(cum_rtn)
  # cum_rtn1 <- cum_rtn
  # for (i in 2:nrow(cum_rtn)){
  #   cum_rtn1[i, ] <- cum_rtn[i, ]-cum_rtn[i-1, ]
  # }
  rtn <- rtn[CT == max(CT) & ZERO_CT <= CT * max_zero_pct]
  rtn <- dcast.data.table(rtn, DATE ~ INNER_CODE, value.var = "CUM_RTN", drop = TRUE)
  rtn <- as.xts.data.table(rtn)
  rtn
}


cal_excess_rtn <- function(cum_rtn, bmk_val) {

  bmk <- merge(cum_rtn[, 0], bmk_val, all = c(TRUE, FALSE))
  bmk <- bmk / matrix(bmk[1, ], nrow = nrow(bmk), ncol(bmk)) - 1
  log(cum_rtn + 1) - log(matrix(bmk, nrow = nrow(cum_rtn), ncol = ncol(cum_rtn)) + 1)
}

opt_ptf <- function(model, candle, date, last_pos) {

  model_piece <- model[EndDate == date, .(INNER_CODE = InnerCode, SIZE, SCORE)]
  setkey(model_piece, INNER_CODE)
  #
  # obj <- setNames(model_piece$SCORE, model_piece$INNER_CODE)
  cum_rtn <- cal_cum_rtn(candle, date, model_piece$INNER_CODE)
  obj <- cum_rtn[nrow(cum_rtn),]/nrow(cum_rtn)


  # cum_rtn <- cal_excess_rtn(cum_rtn, bmk_val)
  # obj <- obj[names(obj) %in% colnames(cum_rtn)]

  stopifnot(is.data.table(last_pos))
  setkey(last_pos, INNER_CODE)
  ini_weight <- last_pos[J(as.integer(colnames(cum_rtn))), GCAMCPUB::na_fill(WEIGHT, 0.0)]
  stopifnot(ini_weight >= 0)
  # obj <- setNames(as.double(cum_rtn[nrow(cum_rtn), ]), colnames(cum_rtn))

  # ptf <- FRAPO::PCDaR(as.matrix(cum_rtn + 1.0), alpha = 0.95, bound = 0.05, softBudget = TRUE)
  # ptf <- ptf@weights
  # ptf <- round(ptf, 6)
  # ptf <- ptf[ptf > 0]
  opt_res <-
    cdar_optimizer_frapo(
      ini_weight,
      obj,
      group,
      cum_rtn = cum_rtn,
      cvar_alpha = 0.95,
      cvar_upper = 0.1,
      group_bound = 0.3,
      max_sec_weight = 0.04,
      max_total_weight = 1,
      # min_total_weight = 0.40
      turnover_bound = 0.3
    )

  ptf <- opt_res$weight
  ptf[, DATE := date]
  ptf[, CDAR := opt_res$cdar]
  ptf
}


opt_ptfs <- function(model, candle, verbose = FALSE) {

  dates <- sort(unique(model$EndDate))
  dates <- dates[dates >= date_from & dates <= date_to]
  # dates <- dates[xts::endpoints(dates, k = 1)]
  ptfs <- vector("list", length(dates))
  prev_pos <- data.table(INNER_CODE = as.integer(model$InnerCode), WEIGHT = as.double(0))
  for (i in seq_along(dates)) {
    date <- dates[i]
    if (i == 1) {
      ptfs[[i]] <- opt_ptf(model, candle, date, prev_pos)
    } else {
      ptfs[[i]] <- opt_ptf(model, candle, date, ptfs[[i - 1]][,1:2])
    }
  }

  # ptfs <- purrr::map(dates, function(date) {
  #   if (isTRUE(verbose)) {
  #     GCAMCPUB::log_info("run cdar optimizer ptf @", date, "...")
  #   }
  #   opt_ptf(model, candle, date)
  # })
  ptfs <- data.table::rbindlist(ptfs)
  setkey(ptfs, DATE, INNER_CODE)
}


log_cumr <- function(ptf, from_to, bmk = c("CSI300", "CSI500", "CSI800")) {

  bmk <- match.arg(bmk)
  bmk_code <- gildata_index_code(bmk)
  from_to <- GCAMCQT::as_from_to(from_to)
  qt <- GCAMCQT::read_gildata("qt_indexquote", from_to = from_to, .load = FALSE)
  bmk_val <- qt$qt_indexquote[J(bmk_code), .(DATE = TRADINGDAY, BMK_VALUE = CLOSEPRICE)]
  setkey(bmk_val, DATE)

  setkey(ptf, DATE, INNER_CODE)
  res <- backtest_by_weight(td, init_capital = 1e10,
                            ptf,
                            from_to = from_to,
                            record_order = FALSE)

  nav <- res$nav[, .(DATE, PTF = LOG_CUMR)]
  setkey(nav, DATE)
  bmk_val_matched <- bmk_val[J(nav$DATE), BMK_VALUE, roll = TRUE]
  nav[, BMK := log(bmk_val_matched) - log(bmk_val_matched[1L])]
  nav[, EXCESS := PTF - BMK]

  as.xts.data.table(nav)
}


rlzd_cdar <- function(dr, alpha) {
  cumr <- cumsum(zoo::na.fill(dr, 0.0))
  dd <- cummax(cumr) - cumr
  mean(dd[dd >= quantile(dd, probs = alpha)])
}



log_cumr2 <- function(ptf, bmk = c("CSI300", "CSI500", "CSI800")) {

  bmk <- match.arg(bmk)
  bmk_code <- gildata_index_code(bmk)
  qt <- GCAMCQT::read_gildata("qt_indexquote", from_to = c(date_from, date_to), .load = FALSE)
  bmk_val <- qt$qt_indexquote[J(bmk_code), .(DATE = TRADINGDAY, BMK_VALUE = CLOSEPRICE)]
  setkey(bmk_val, DATE)

  setkey(ptf, DATE, INNER_CODE)
  res <- backtest_by_weight(td, init_capital = 1e10,
                            ptf[DATE >= GCAMCPUB::f_date_begin(date_from, "month") & DATE <= date_to],
                            from_to = c(date_from, date_to),
                            record_order = FALSE)

  nav <- res$nav[, .(DATE, PTF = LOG_CUMR)]
  setkey(nav, DATE)
  bmk_val_matched <- bmk_val[J(nav$DATE), BMK_VALUE, roll = TRUE]
  nav[, BMK := log(bmk_val_matched) - log(bmk_val_matched[1L])]
  nav[, EXCESS := PTF - BMK]

  as.xts.data.table(nav)
}


reblance <- function(dat) {

  dat - matrix(dat[1, ], nrow = nrow(dat), ncol = ncol(dat), byrow = TRUE)
}


bmk_ptf <- function(bmk) {

  bmk_size_dat <- bmk_size(bmk, model)

  size_constr <- function(date) {
    bmk_size_dat[["size"]][J(date), SIZE] * c(0.95, 1.05)
  }

  stock_wt <- function(date, default_max_wt) {

    model_inner_code <- model[EndDate == date, InnerCode]

    bmk_wt <- bmk_size_dat$raw[DATE == date, .(INNER_CODE, WEIGHT)]
    setkey(bmk_wt, INNER_CODE)
    univ_wt <- bmk_wt[J(model_inner_code), WEIGHT]
    upper <- list(
      ind = 1:length(model_inner_code),
      val = local({
        tmp <- univ_wt
        pos <- which(!is.na(tmp))
        tmp[pos] <- pmax(pmin(tmp[pos] * 10,  0.3), default_max_wt)
        tmp[-pos] <- default_max_wt
        tmp
      })
    )

    list(upper = upper)
  }


  dates <- sort(unique(model$EndDate))
  ptfs <- purrr::map(dates, ~size_controlled_ptf(model, ., 0.05, size_constr(.), bounds = stock_wt(., 0.05)))
  ptfs <- data.table::rbindlist(ptfs)
  ptfs
}


size_ptf <- function(size_range, from_to) {

  size_constr <- function(date) {
    size_range
  }

  dates <- sort(unique(model$EndDate))
  dates <- dates[dates >= from_to[1L] & dates <= from_to[2L]]
  ptfs <- purrr::map(dates, ~size_controlled_ptf(model, ., 0.05, size_constr(.)))
  ptfs <- data.table::rbindlist(ptfs)
  log_cumr2(ptfs, "CSI300")
}
