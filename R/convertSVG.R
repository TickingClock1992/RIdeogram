##' convert svg to png or other format
##'
##'
##' @title convertSVG
##' @param svg svg file
##' @param device target format
##' @param width output width
##' @param height output height
##' @param dpi output dpi
##' @return invisible grob object
##' @importFrom rsvg rsvg_svg
##' @importFrom grImport2 readPicture
##' @importFrom grImport2 pictureGrob
##' @importFrom ggplot2 ggsave
##' @export
##' @author Zhaodong Hao, Dekang Lv, Guangchuang Yu, Ying Ge, Jisen Shi, Jinhui Chen
convertSVG <- function(svg, device = "png", width = 8.2677, height = 11.6929, dpi = 300) {
  f <- tempfile(fileext = ".svg")
  rsvg::rsvg_svg(svg, f)
  x <- grImport2::readPicture(f)
  g <- grImport2::pictureGrob(x)
  path <- getwd()
  ggsave(g, filename = paste("chromosome.", device, sep = ""), device = device, path = path, width = width, height = height, dpi = dpi)
  invisible(g)
}
