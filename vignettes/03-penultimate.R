## ----penult, fig.width = 10, out.width = '90%', fig.align = "center"----------
library(mev)
set.seed(123)
x  <- seq(1, 30, length = 200)
fitted <- t(replicate(n = 1000,
                      fit.gev(apply(
                        matrix(rlnorm(30 * 100), ncol = 30), 1, max
                      ),
                      method = "BFGS")$est))
penult.bm.lnorm <-
  mev::smith.penult(model = "bm",
                    m = (m <- 30),
                    family = "lnorm")

par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
hist(
  fitted[, 3],
  probability = TRUE,
  breaks = 20,
  xlab = "estimated shape",
  main = ""
)
segments(
  y0 = -0.5,
  y1 = 0,
  x0 = penult.bm.lnorm$shape,
  col = "red",
  lwd = 2
)
segments(
  y0 = -0.5,
  y1 = 0,
  x0 = 0,
  col = "black",
  lwd = 2
)

p30 <- penult.bm.lnorm
x = seq(0, 100, length = 400)
N <- 1000
N30 = N / 30
plot(
  x,
  N * exp((N - 1) * plnorm(x, log.p = TRUE)) *
    dlnorm(x),
  type = "l",
  bty = "l",
  ylim = c(0, 0.1),
  xlab = "x",
  ylab = "density"
)
# Get parameters of maximum of N GEV
maxstabp <- function(loc, scale, shape, N) {
  if (!isTRUE(all.equal(shape, 0, check.attributes = FALSE))) {
    mut <- loc - scale * (1 - N ^ shape) / shape
    sigmat = scale * N ^ shape
    return(c(
      loc = mut,
      scale = sigmat,
      shape = shape
    ))
  } else{
    mut <- loc + scale * log(N)
    return(c(
      loc = mut,
      scale = scale,
      shape = shape
    ))
  }
}
p30e <-
  maxstabp(
    loc = p30$loc,
    scale = p30$scale,
    shape = p30$shape,
    N = N30
  )
lines(
  x,
  evd::dgev(
    x,
    loc = p30e['loc'],
    scale = p30e['scale'],
    shape = p30e['shape']
  ),
  col = "red",
  lwd = 2,
  lty = 2
)
p30l <-
  maxstabp(
    loc = p30$loc,
    scale = p30$scale,
    shape = 0,
    N = N30
  )
lines(
  x,
  evd::dgev(
    x,
    loc = p30e['loc'],
    scale = p30l['scale'],
    shape = 0
  ),
  col = 1,
  lty = 3,
  lwd = 1
)
legend(
  x = "topright",
  legend = c("exact", "ultimate", "penultimate"),
  col = c(1, 1, 2),
  lwd = c(2, 2, 1),
  bty = "n",
  lty = c(1, 3, 2)
)

## ----penultimateshape, fig.width = 10, fig.height = 10, fig.align = "center", out.width = '70%'----
u <- seq(0.8, 0.9999, by = 0.0001)
penult <- smith.penult(family = "norm",
                       method = "pot",
                       u = qnorm(u))
plot(
  u,
  penult$shape,
  bty = "l",
  type = 'l',
  xlab = 'quantile',
  ylab = 'penultimate shape'
)

