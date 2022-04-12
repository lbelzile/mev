## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ----tstabplots, out.width = '90%', fig.width = 8, fig.height = 10, fig.align = "center"----
library(mev)
data(maiquetia, package = "mev")
day <- seq.Date(from = as.Date("1961-01-01"), 
                to = as.Date("1999-12-31"), 
                by = "day")
nzrain <- maiquetia[substr(day, 1, 4) < 1999 & maiquetia > 0]
# Keep non-zero rainfall, exclude 1999 observations
nzrain <- maiquetia[substr(day, 3, 4) < 99 & maiquetia > 0]
thcan <- quantile(nzrain, 
                  seq(0.8, 0.99, length.out = 25))
tstab.gpd(xdat = nzrain, 
          thresh = thcan, 
          method = "profile")

## ----NCdiag, out.width = '90%',fig.width = 10, fig.height = 6, fig.align = "center", out.width = '90%'----
fNCdiag <- NC.diag(x = nzrain, u = thcan)

## ----tstagegp, out.width = '90%', fig.width = 8, fig.height = 15, fig.align = "center"----
tstab.egp(xdat = nzrain, 
          thresh = thcan, 
          model = "egp2")

## ----Wdiag, out.width = '90%',fig.width = 10, fig.height = 12, fig.align = "center"----
fWdiag <-
  W.diag(
    xdat = nzrain,
    model = "nhpp",
    u = thcan,
    plots = c("WN", "PS")
  )

