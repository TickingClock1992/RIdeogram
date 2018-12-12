##' convert svg to pdf
##'
##'
##' @title svg2pdf
##' @param svg svg file
##' @param pdf pdf file
##' @param width output width
##' @param height output height
##' @return invisible grob object
##' @importFrom rsvg rsvg_svg
##' @importFrom grImport2 readPicture
##' @importFrom grImport2 pictureGrob
##' @importFrom ggplot2 ggsave
##' @export
##' @author Guangchuang Yu
svg2pdf <- function(svg, pdf, width = 210, height = 297) {
  f <- tempfile(fileext = ".svg")
  rsvg::rsvg_svg(svg, f)
  x <- grImport2::readPicture(f)
  g <- grImport2::pictureGrob(x)
  path <- getwd()
  ggsave(g, filename = "chromosome.pdf", path = path, width = width, height = height, units = "mm")
  invisible(g)
}
