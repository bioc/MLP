#' Plot Quantile Curves for an MLP Run 
#' @param x0 TODO
#' @param y0 TODO
#' @param hqi TODO
#' @param xtp TODO
#' @param qi TODO
#' @param lqi TODO
#' @param sym TODO defaults to TRUE
#' @param main main title for the quantile curves plot; defaults to NULL; if NULL
#' no title is added
#' @return no return value; a quantile curve plot is drawn to the current device
#' @examples pathExampleMLPResult <- system.file("exampleFiles", "exampleMLPResult.rda", package = "MLP")
#' load(pathExampleMLPResult)
#' @noRd
plotQuantileCurves <- function(x0, y0, hqi, xtp, qi, lqi, sym = TRUE, main = NULL) {
 
  mainTitle <- if (is.null(main)) "" else main
  
  darkBlue <- "#08306B"
  someBlue <- "#2171B5"
  someLightBlue <- "#539ECC"

  plot(x0, y0, xlab = "n", ylab = "MLP", axes = FALSE, col = darkBlue, pch = ".",
      main = mainTitle)
  axis(2, lwd = 1.5, col = darkBlue)
  atPositions <- axis(1, labels = FALSE)
  axis(1, lwd = 1.5, at = atPositions, labels = atPositions^2, las = 1, col = darkBlue)
  op <- par(cex = 0.5)
  npc <- 16
  
  if (!is.null(hqi)){ 
    i2 <- y0 > xtp[,length(qi)]
    points(x0[i2], y0[i2], col = someBlue, pch = npc) 
    if (length(i2[i2 == TRUE]) > 0){
      text(x = x0[i2], y = jitter(y0[i2], factor=4), labels = names(x0[i2]), cex = 1.25, col = darkBlue)
    }
  }
  if (!is.null(lqi)){ 
    i2 <- y0 < xtp[,1] 
    points(x0[i2], y0[i2], col = someBlue, pch = npc) 
  }
  if (sym){ 
    i2 <- y0 < -xtp[, length(qi)] 
    points(x0[i2], y0[i2], col = someBlue,pch = npc) 
  }
  par(op)
  
  matlines(x0[i2 <- sort.list(x0)], xtp[i2,], col = someLightBlue, lty = "solid")
  if (sym) 
    matlines(x0[i2], -xtp[i2,], col = 3, lty = "solid")
}
