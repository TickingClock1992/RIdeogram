##' convert svg to png or other format
##'
##'
##' @title convertSVG
##' @param svg svg file
##' @param file output file name
##' @param device target format
##' @param width output width
##' @param height output height
##' @param dpi output dpi
##' @return invisible grob object
##' @importFrom rsvg rsvg_svg
##' @importFrom grImport2 readPicture
##' @importFrom grImport2 pictureGrob
##' @importFrom ggplot2 ggsave
##' @importFrom tools file_ext
##' @export
##' @rdname convertSVG
##' @author Zhaodong Hao, Dekang Lv, Guangchuang Yu, Ying Ge, Jisen Shi, Jinhui Chen
convertSVG <- function(svg, file = "chromosome", device = NULL, width = 8.2677, height = 11.6929, dpi = 300) {
  f <- tempfile(fileext = ".svg")
  rsvg::rsvg_svg(svg, f)
  x <- grImport2::readPicture(f)
  g <- grImport2::pictureGrob(x)
  path <- getwd()
  ext <- file_ext(file)
  outfile <- file
  if (is.null(device)) {
      if (ext == "") {
          stop("please specifiy 'device' to one of 'pdf', 'tiff', 'png', or 'jpg'.")
      } else {
          device <- ext
      }
  }

  if (! tolower(device) %in% c('pdf', 'tiff', 'png', 'jpg')) {
      stop("'device' should be one of 'pdf', 'tiff', 'png', or 'jpg'.")
  }

  if (file_ext(file) != device)
      outfile <- paste0(outfile, ".", device)

  ggsave(g, filename = outfile, device = device, path = path, width = width, height = height, dpi = dpi)
  invisible(g)
}

##' @rdname convertSVG
##' @export
svg2pdf <- function(svg, file, width = 8.2677, height = 11.6929, dpi = 300) {
    convertSVG(svg, file, device = 'pdf', width = width, height = height, dpi = dpi)
}

##' @rdname convertSVG
##' @export
svg2png <- function(svg, file, width = 8.2677, height = 11.6929, dpi = 300) {
    convertSVG(svg, file, device = 'png', width = width, height = height, dpi = dpi)
}

##' @rdname convertSVG
##' @export
svg2tiff <- function(svg, file, width = 8.2677, height = 11.6929, dpi = 300) {
    convertSVG(svg, file, device = 'tiff', width = width, height = height, dpi = dpi)
}

##' @rdname convertSVG
##' @export
svg2jpg <- function(svg, file, width = 8.2677, height = 11.6929, dpi = 300) {
    convertSVG(svg, file, device = 'jpg', width = width, height = height, dpi = dpi)
}

