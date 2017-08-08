library(Rglpk)
library(dplyr)
library(tidyr)
library(dtplyr)
library(data.table)
library(lubridate)
library(GCAMCBT)
source("cvar_public_functions.R")

date_from <- ymd(20091231)
date_to <- ymd(20170717)

td <- trading_desk()
candle <- candle_ashare(c(date_from, date_to))
load_candle(td, candle)

# read model data
model <- GCAMCPUB::readDtRds("model_20170802.rds")
group <- readr::read_rds("E:/RWD/Deng/cdar/group.rds")
invisible(add_size_factor(model))
invisible(normalize_model_score(model))

stop()
ptfs <- opt_ptfs(model[EndDate >= date_from & SCORE > 1], candle, verbose = TRUE)

ptfs <- local({
  dates <- sort(unique(model$EndDate))
  dates <- dates[xts::endpoints(dates)]
  opt_ptfs(model[SCORE > 1 & EndDate %in% dates], candle, verbose = TRUE)
})

res <- log_cumr(ptfs, c(date_from, date_to), "CSI300")

local({

  dr <- xts::diff.xts(res[, "PTF"])
  cdar095 <- function(x) rlzd_cdar(x, 0.95)
  cdar_res <- zoo::rollapplyr(dr, width = 90L, FUN = cdar095)
  colnames(cdar_res) <- "CDaR095"
  purrr::walk(c(0.05, 0.10, 0.15), ~print(sum(cdar_res >= ., na.rm = TRUE) / sum(!is.na(cdar_res))))
  plot(cbind(res, cdar_res, 0.05, 0.10, 0.15))
})

rlzd_cdar(xts::diff.xts(res["201607/201609", "PTF"]), alpha = 0.95)

plot(reblance(res["2006/2017"]), legend.loc = "topleft")
plot(reblance(res["201411/201603"]), legend.loc = "topleft")
plot(reblance(res["20160105/201612"]), legend.loc = "topleft")
plot(reblance(res["20160105/201707"]), legend.loc = "topleft")
plot(reblance(res["2017/201707"]), legend.loc = "topleft")
plot(reblance(res["20170306/201707"]), legend.loc = "topleft")
ptfs[, .(WEIGHT = sum(WEIGHT)), keyby = DATE] %>% as.xts.data.table() %>% plot(.)
model[ptfs[DATE == max(DATE)], on = c("InnerCode" = "INNER_CODE", "EndDate" = "DATE")]

local({

  cumr <- reblance(res)[, c("PTF", "BMK")]
  dr <- xts::diff.xts(cumr)
  PerformanceAnalytics::charts.PerformanceSummary(dr)
  PerformanceAnalytics::InformationRatio(dr[, "PTF"], dr[, "BMK"]) %>% print(.)
  PerformanceAnalytics::SharpeRatio.annualized(dr[, "PTF"]) %>% print(.)
})


local({

  model[ptfs, on = c("InnerCode" = "INNER_CODE", "EndDate" = "DATE")][
    , .(SIZE = weighted.mean(SIZE, WEIGHT, na.rm = TRUE)), keyby = EndDate] %>%
    as.xts.data.table() %>%
    plot(.)
})

# dt <- data.table(a = 1:26, b = LETTERS, value = 1)
# dcast.data.table(dt, b ~ a, value.var = "value", fill = 0)
