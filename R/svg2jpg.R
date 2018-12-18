##' convert svg to jpg
##'
##'
##' @title svg2jpg
##' @param svg svg file
##' @param jpg jpg file
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
svg2jpg <- function(svg, jpg, width = 210, height = 297, dpi = 300) {
  f <- tempfile(fileext = ".svg")
  rsvg::rsvg_svg(svg, f)
  x <- grImport2::readPicture(f)
  g <- grImport2::pictureGrob(x)
  path <- getwd()
  ggsave(g, filename = "chromosome.jpg", device = "jpeg", path = path, width = width, height = height, units = "mm", dpi = dpi)
  invisible(g)
}
