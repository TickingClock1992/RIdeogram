##' ideogram with overlaid heatmap annotation and optional track label
##'
##'
##' @title ideogram
##' @param karyotype karyotype data
##' @param overlaid overlaid data
##' @param label track label data
##' @param label_type track label type, only support four types: marker, heatmap, line and polygon
##' @param synteny synteny data
##' @param colorset1 overlaid heatmap color
##' @param colorset2 label heatmap color
##' @param width width of plot region
##' @param Lx position of legend (x)
##' @param Ly position of legend (y)
##' @param output output file, only svg is supported
##' @return output file
##' @importFrom grDevices colorRampPalette
##' @importFrom scales rescale
##' @export
##' @examples
##' # Loading the package
##' require(RIdeogram)
##'
##' # Loading the testing data
##' data(human_karyotype, package="RIdeogram")
##' data(gene_density, package="RIdeogram")
##' data(Random_RNAs_500, package="RIdeogram")
##'
##' # Checking the data format
##' head(human_karyotype)
##' head(gene_density)
##' head(Random_RNAs_500)
##'
##' # Running the function
##' ideogram(karyotype = human_karyotype)
##' convertSVG("chromosome.svg", device = "png")
##'
##' # Then, you will find a SVG file and a PNG file in your Working Directory.
##'
##' @author Zhaodong Hao, Dekang Lv, Ying Ge, Jisen Shi, Dolf Weijers, Guangchuang Yu, Jinhui Chen
##'
ideogram <- function(karyotype, overlaid = NULL, label = NULL, label_type = NULL, synteny = NULL, colorset1 = c("#4575b4", "#ffffbf", "#d73027"), colorset2 = c("#b35806", "#f7f7f7", "#542788"), width = 170, Lx = 160, Ly = 35, output = "chromosome.svg") {
  karyotype <- karyotype
  mydata <- overlaid
  mydata_interval <- label
  synteny_data <- synteny

  col_num <- ncol(karyotype)

  if (!is.null(synteny)) {
    if (length(table(karyotype$species)) == 2){
      # karyotype rect
      # 1.5cm was left for the corresponding species name & interval between chr is 0.1cm
      chr_width_total_1 <- ((17 - 1.5) - (table(karyotype$species)[1] - 1) * 0.1) * 35.43307
      karyotype_total_1 <- sum(karyotype$End[1:table(karyotype$species)[1]])

      for (i in 1:table(karyotype$species)[1]) {
        j = i - 1
        if (i == 1)
          karyotype[i, 8] <- 3.5*35.43307
        else
          karyotype[i, 8] <- 3.5*35.43307 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + (i-1)*0.1*35.43307
      }
      names(karyotype)[8] <- "x"

      karyotype$y = 2.5*35.43307

      # karyotype rect

      # interval between chr is 0.1cm
      chr_width_total_2 <- ((17 - 1.5) - (table(karyotype$species)[2] - 1) * 0.1) * 35.43307
      karyotype_total_2 <- sum(karyotype$End[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])])

      for (i in (table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])) {
        karyotype[i, 9] = (2.5 + 0.5 + 5)*35.43307
        j = i - 1
        if (i == table(karyotype$species)[1] + 1)
          karyotype[i, 8] <- 3.5*35.43307
        else
          karyotype[i, 8] <- 3.5*35.43307 + sum(karyotype[(table(karyotype$species)[1] + 1):j,3]) * chr_width_total_2 / karyotype_total_2 + (i-table(karyotype$species)[1]-1)*0.1*35.43307
      }

      #
      karyotype$rx = 2

      karyotype$ry = 2

      for (i in 1:table(karyotype$species)[1]) {
        karyotype[i, 12] <- karyotype[i, 3] * chr_width_total_1 / karyotype_total_1
      }
      for (i in (table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])) {
        karyotype[i, 12] <- karyotype[i, 3] * chr_width_total_2 / karyotype_total_2
      }
      names(karyotype)[12] <- "width"

      karyotype$rect_outer = paste("<rect x=\"", karyotype$x, "\" y=\"", karyotype$y, "\" rx=\"", karyotype$rx,  "\" ry=\"", karyotype$ry,  "\" width=\"", karyotype$width,  "\" height=\"17.71654\" style=\"fill:#f7f7f7;stroke:black;stroke-width:0.5\"/>", sep = "")
      karyotype$rect_inner = paste("<rect x=\"", karyotype$x, "\" y=\"", karyotype$y + 3, "\" rx=\"", karyotype$rx,  "\" ry=\"", karyotype$ry,  "\" width=\"", karyotype$width,  "\" height=\"11.71654\" style=\"fill:#", karyotype$fill, ";stroke:#", karyotype$fill, ";stroke-width:0.5\"/>", sep = "")

      #synteny
      for (i in 1: nrow(synteny)) {
        j = synteny[i, 1] - 1
        if (synteny[i, 1] == "1")
          synteny[i, 8] <- 3.5*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1
        else
          synteny[i, 8] <- 3.5*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.1*35.43307
      }
      names(synteny)[8] <- "x1"

      synteny$y1 <- (2.5 + 0.5)*35.43307

      for (i in 1: nrow(synteny)) {
        j = synteny[i, 4] - 1
        if (synteny[i, 4] == "1")
          synteny[i, 10] <- 3.5*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2
        else
          synteny[i, 10] <- 3.5*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.1*35.43307
      }
      names(synteny)[10] <- "x2"

      synteny$y2 <- (2.5 + 0.5 + 5)*35.43307

      for (i in 1: nrow(synteny)) {
        j = synteny[i, 4] - 1
        if (synteny[i, 4] == "1")
          synteny[i, 12] <- 3.5*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2
        else
          synteny[i, 12] <- 3.5*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.1*35.43307
      }
      names(synteny)[12] <- "x3"

      synteny$y3 <- (2.5 + 0.5 + 5)*35.43307

      for (i in 1: nrow(synteny)) {
        j = synteny[i, 1] - 1
        if (synteny[i, 1] == "1")
          synteny[i, 14] <- 3.5*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1
        else
          synteny[i, 14] <- 3.5*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.1*35.43307
      }
      names(synteny)[14] <- "x4"

      synteny$y4 <- (2.5 + 0.5)*35.43307

      synteny$x5 <- synteny$x1

      synteny$y5 <- synteny$y1 + 1.5*35.43307

      synteny$x6 <- (synteny$x1 + synteny$x2) / 2

      synteny$y6 <- (synteny$y1 + synteny$y2) / 2

      synteny$x7 <- synteny$x3

      synteny$y7 <- synteny$y3 - 1.5*35.43307

      synteny$x8 <- (synteny$x3 + synteny$x4) / 2

      synteny$y8 <- (synteny$y3 + synteny$y4) / 2

      synteny$path <- paste(paste("<path d=\"M", synteny$x1, synteny$y1, "Q", synteny$x5, synteny$y5, synteny$x6, synteny$y6, "T", synteny$x2, synteny$y2, "L", synteny$x3, synteny$y3, "Q", synteny$x7, synteny$y7, synteny$x8, synteny$y8, "T", synteny$x4, synteny$y4, "Z\""), paste(" stroke=\"#", synteny$fill, "\" fill=\"#", synteny$fill, "\" style=\"stroke-width: 0.5px;\"/>", sep = ""), sep = "")

      # text model <text x y font-size fill >words</text>

      species <- data.frame(paste("<text x=\"", 2 * 35.43307, "\" y=\"", (2.5 + 0.5) * 35.43307 - 3, "\" font-size=\"", karyotype[1,6], "\" fill=\"#", karyotype[1,7], "\" >", karyotype[1,5], "</text>", sep = ""),
                            paste("<text x=\"", 2 * 35.43307, "\" y=\"", (2.5 + 0.5 + 5 + 0.5) * 35.43307 - 3, "\" font-size=\"", karyotype[table(karyotype$species)[1]+1,6], "\" fill=\"#", karyotype[table(karyotype$species)[1]+1,7], "\" >", karyotype[table(karyotype$species)[1]+1,5], "</text>", sep = ""),
                            stringsAsFactors = FALSE
      )

      for (i in 1 : table(karyotype$species)[1]){
        karyotype[i,15] = paste("<text x=\"", (karyotype[i, 8] + karyotype[i, 8] + karyotype[i, 12])/2 - nchar(karyotype[i, 1]) * 2.1, "\" y=\"", karyotype[i, 9] - 0.1 * 35.43307, "\" font-size=\"9\" fill=\"black\" >", karyotype[i, 1], "</text>", sep = "")
      }
      for (i in (table(karyotype$species)[1] + 1) : (table(karyotype$species)[1] + table(karyotype$species)[2])){
        karyotype[i,15] = paste("<text x=\"", (karyotype[i, 8] + karyotype[i, 8] + karyotype[i, 12])/2 - nchar(karyotype[i, 1]) * 2.1, "\" y=\"", karyotype[i, 9] + (0.1 + 0.5 + 0.15) * 35.43307, "\" font-size=\"9\" fill=\"black\" >", karyotype[i, 1], "</text>", sep = "")
      }
      names(karyotype)[15] <- "text"

    } else if (length(table(karyotype$species)) == 3){
      #karyotype_1
      chr_width_total_1 <- (11 - (table(karyotype$species)[1] - 1) * 0.04) * 35.43307
      karyotype_total_1 <- sum(karyotype$End[1:table(karyotype$species)[1]])

      for (i in 1:table(karyotype$species)[1]) {
        j = i - 1
        if (i == 1)
          karyotype[i, 8] <- 3*35.43307
        else
          karyotype[i, 8] <- 3*35.43307 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + (i-1)*0.04*35.43307
      }
      names(karyotype)[8] <- "x"

      karyotype$y = (2.5 + 11.25833)*35.43307

      karyotype$rx = 1

      karyotype$ry = 1

      for (i in 1:table(karyotype$species)[1]) {
        karyotype[i, 12] <- karyotype[i, 3] * chr_width_total_1 / karyotype_total_1
      }
      names(karyotype)[12] <- "width"

      for (i in 1:table(karyotype$species)[1]) {
        karyotype[i, 13] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9], "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"8\" style=\"fill:f7f7f7;stroke:black;stroke-width:0.1;fill-opacity:0\"/>", sep = "")
        karyotype[i, 14] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9] + 2, "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"4\" style=\"fill:#", karyotype[i, 4], ";stroke:#", karyotype[i, 4], ";stroke-width:0.1\"/>", sep = "")
      }
      names(karyotype)[13] <- "rect_outer"
      names(karyotype)[14] <- "rect_inner"

      # karyotype_2
      chr_width_total_2 <- (11 - (table(karyotype$species)[2] - 1) * 0.04) * 35.43307
      karyotype_total_2 <- sum(karyotype$End[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])])

      for (i in (table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])) {
        j = i - 1
        if (i == table(karyotype$species)[1] + 1)
          karyotype[i, 8] <- 3*35.43307
        else
          karyotype[i, 8] <- 3*35.43307 + sum(karyotype[(table(karyotype$species)[1] + 1):j,3]) * chr_width_total_2 / karyotype_total_2 + (i- table(karyotype$species)[1] - 1)*0.04*35.43307
      }
      names(karyotype)[8] <- "x"

      karyotype$y = (2.5 + 11.25833)*35.43307

      for (i in (table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])) {
        karyotype[i, 12] <- karyotype[i, 3] * chr_width_total_2 / karyotype_total_2
      }
      names(karyotype)[12] <- "width"

      for (i in (table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2])) {
        karyotype[i, 13] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9], "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"8\" transform=\"rotate(-60, 70.86614 487.4999)\" style=\"fill:#f7f7f7;stroke:black;stroke-width:0.1;fill-opacity:0\"/>", sep = "")
        karyotype[i, 14] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9] + 2, "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"4\" transform=\"rotate(-60, 70.86614 487.4999)\" style=\"fill:#", karyotype[i, 4], ";stroke:#", karyotype[i, 4], ";stroke-width:0.1\"/>", sep = "")
      }

      # karyotype_3
      chr_width_total_3 <- (11 - (table(karyotype$species)[3] - 1) * 0.04) * 35.43307
      karyotype_total_3 <- sum(karyotype$End[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + table(karyotype$species)[3])])

      for (i in (table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + table(karyotype$species)[3])) {
        j = i -1
        if (i == table(karyotype$species)[1] + table(karyotype$species)[2] + 1)
          karyotype[i, 8] <- 3*35.43307
        else
          karyotype[i, 8] <- 3*35.43307 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):j,3]) * chr_width_total_3 / karyotype_total_3 + (i- table(karyotype$species)[1] - table(karyotype$species)[2] - 1)*0.04*35.43307
      }
      names(karyotype)[8] <- "x"

      karyotype$y = (2.5 + 11.25833)*35.43307

      for (i in (table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + table(karyotype$species)[3])) {
        karyotype[i, 12] <- karyotype[i, 3] * chr_width_total_3 / karyotype_total_3
      }
      names(karyotype)[12] <- "width"

      for (i in (table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + table(karyotype$species)[3])) {
        karyotype[i, 13] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9], "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"8\" transform=\"rotate(60, 531.496 487.4999)\" style=\"fill:#f7f7f7;stroke:black;stroke-width:0.1;fill-opacity:0\"/>", sep = "")
        karyotype[i, 14] <- paste("<rect x=\"", karyotype[i, 8], "\" y=\"", karyotype[i, 9] + 2, "\" rx=\"", karyotype[i, 10],  "\" ry=\"", karyotype[i, 11],  "\" width=\"", karyotype[i, 12],  "\" height=\"4\" transform=\"rotate(60, 531.496 487.4999)\" style=\"fill:#", karyotype[i, 4], ";stroke:#", karyotype[i, 4], ";stroke-width:0.1\"/>", sep = "")
      }

      #synteny

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 9] <- 3*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1
          else
            synteny[i, 9] <- 3*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.04*35.4330
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 9] <- 3*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1
          else
            synteny[i, 9] <- 3*35.43307 + synteny[i, 2] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.04*35.4330
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 9] <- (1*35.43307 + synteny[i, 2] * chr_width_total_2 / karyotype_total_2) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
          else
            synteny[i, 9] <- (1*35.43307 + synteny[i, 2] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
        }
      }
      names(synteny)[9] <- "x1"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          synteny[i, 10] <- (2.5 + 11.25833)*35.43307
        } else if (synteny[i, 8] == 2){
          synteny[i, 10] <- (2.5 + 11.25833)*35.43307
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 10] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 2] * chr_width_total_2 / karyotype_total_2) * sin(pi*(1/3)) + 4
          else
            synteny[i, 10] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 2] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * sin(pi*(1/3)) + 4
        }
      }
      names(synteny)[10] <- "y1"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 11] <- 3*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1
          else
            synteny[i, 11] <- 3*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.04*35.4330
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 11] <- 3*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1
          else
            synteny[i, 11] <- 3*35.43307 + synteny[i, 3] * chr_width_total_1 / karyotype_total_1 + sum(karyotype[1:j,3]) * chr_width_total_1 / karyotype_total_1 + j*0.04*35.4330
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 11] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
          else
            synteny[i, 11] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
        }
      }
      names(synteny)[11] <- "x2"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          synteny[i, 12] <- (2.5 + 11.25833)*35.43307
        } else if (synteny[i, 8] == 2){
          synteny[i, 12] <- (2.5 + 11.25833)*35.43307
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 12] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3)) * sin(pi*(1/3)) + 4
          else
            synteny[i, 12] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * sin(pi*(1/3)) + 4
        }
      }
      names(synteny)[12] <- "y2"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 13] <- (1*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
          else
            synteny[i, 13] <- (1*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 13] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
          else
            synteny[i, 13] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 13] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
          else
            synteny[i, 13] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
        }
      }
      names(synteny)[13] <- "x3"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2) * sin(pi*(1/3)) + 4
          else
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * sin(pi*(1/3)) + 4
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3)) * sin(pi*(1/3)) + 4
          else
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * sin(pi*(1/3)) + 4
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3)) * sin(pi*(1/3)) + 4
          else
            synteny[i, 14] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 6] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * sin(pi*(1/3)) + 4
        }
      }
      names(synteny)[14] <- "y3"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 4] - 1
          if (synteny[i,4] == "1")
            synteny[i, 15] <- (1*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
          else
            synteny[i, 15] <- (1*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 15] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
          else
            synteny[i, 15] <- 15 * 35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * cos(pi*(1/3)) - 8 * cos(pi*(1/6))
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 15] <- (1*35.43307 + synteny[i, 3] * chr_width_total_2 / karyotype_total_2) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
          else
            synteny[i, 15] <- (1*35.43307 + synteny[i, 3] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * cos(pi*(1/3)) + 2 * 35.43307 + 8 * cos(pi*(1/6))
        }
      }
      names(synteny)[15] <- "x4"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2) * sin(pi*(1/3)) + 4
          else
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * sin(pi*(1/3)) + 4
        } else if (synteny[i, 8] == 2){
          j = synteny[i, 4] - 1
          if (synteny[i, 4] == "1")
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3)) * sin(pi*(1/3)) + 4
          else
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (13 * 35.43307 - (1*35.43307 + synteny[i, 5] * chr_width_total_3 / karyotype_total_3 + sum(karyotype[(table(karyotype$species)[1] + table(karyotype$species)[2] + 1):(table(karyotype$species)[1] + table(karyotype$species)[2] + j),3]) * chr_width_total_3 / karyotype_total_3 + j*0.04*35.43307)) * sin(pi*(1/3)) + 4
        } else if (synteny[i, 8] == 3){
          j = synteny[i, 1] - 1
          if (synteny[i, 1] == "1")
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 3] * chr_width_total_2 / karyotype_total_2) * sin(pi*(1/3)) + 4
          else
            synteny[i, 16] <- (2.5 + 11.25833)*35.43307 - (1*35.43307 + synteny[i, 3] * chr_width_total_2 / karyotype_total_2 + sum(karyotype[(table(karyotype$species)[1] + 1):(table(karyotype$species)[1] + j),3]) * chr_width_total_2 / karyotype_total_2 + j*0.04*35.43307) * sin(pi*(1/3)) + 4
        }
      }
      names(synteny)[16] <- "y4"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          synteny[i, 17] <- (synteny[i, 9] + 11 * 35.43307) / 2 - 30
          synteny[i, 18] <- (2.5 + 11.25833)*35.43307 - (synteny[i, 17] - 2 * 35.43307) * tan(pi*(1/6))
          synteny[i, 19] <- (synteny[i, 11] + 11 * 35.43307) / 2 - 30
          synteny[i, 20] <- (2.5 + 11.25833)*35.43307 - (synteny[i, 19] - 2 * 35.43307) * tan(pi*(1/6))
        } else if (synteny[i, 8] == 2){
          synteny[i, 17] <- (synteny[i, 9] + 11 * 35.43307) / 2 + 30
          synteny[i, 18] <- (2.5 + 11.25833)*35.43307 - (15 * 35.43307 - synteny[i, 17]) * tan(pi*(1/6))
          synteny[i, 19] <- (synteny[i, 11] + 11 * 35.43307) / 2 + 30
          synteny[i, 20] <- (2.5 + 11.25833)*35.43307 - (15 * 35.43307 - synteny[i, 19]) * tan(pi*(1/6))
        } else if (synteny[i, 8] == 3){
          synteny[i, 17] <- 8.5 * 35.43307
          synteny[i, 18] <- (synteny[i, 10] + (2.5 + 11.25833 / 2)*35.43307 ) / 2 - 30
          synteny[i, 19] <- 8.5 * 35.43307
          synteny[i, 20] <- (synteny[i, 10] + (2.5 + 11.25833 / 2)*35.43307 ) / 2 - 30
        }
      }
      names(synteny)[17] <- "x5"
      names(synteny)[18] <- "y5"
      names(synteny)[19] <- "x6"
      names(synteny)[20] <- "y6"

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 15], synteny[i, 16], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 9], synteny[i, 10], "L", synteny[i, 11], synteny[i, 12], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 13], synteny[i, 14], "Z\""), paste(" stroke=\"#", synteny[i, 7], "\" fill=\"#", synteny[i, 7], "\" style=\"stroke-width: 0.1px;\"/>", sep = ""), sep = "")
        } else if (synteny[i, 8] == 2){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 15], synteny[i, 16], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 9], synteny[i, 10], "L", synteny[i, 11], synteny[i, 12], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 13], synteny[i, 14], "Z\""), paste(" stroke=\"#", synteny[i, 7], "\" fill=\"#", synteny[i, 7], "\" style=\"stroke-width: 0.1px;\"/>", sep = ""), sep = "")
        } else if (synteny[i, 8] == 3){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 9], synteny[i, 10], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 11], synteny[i, 12], "L", synteny[i, 13], synteny[i, 14], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 15], synteny[i, 16], "Z\""), paste(" stroke=\"#", synteny[i, 7], "\" fill=\"#", synteny[i, 7], "\" style=\"stroke-width: 0.1px;\"/>", sep = ""), sep = "")
        }
      }
      names(synteny)[21] <- "path"


      species <- data.frame(paste("<text x=\"", 8.5 * 35.43307 - nchar(names(table(karyotype[,5]))[1]) * 2.1, "\" y=\"", (2.5 + 11.25833 + 0.8) * 35.43307, "\" font-size=\"", karyotype[1,6], "\" fill=\"#", karyotype[1,7], "\" >", names(table(karyotype[,5]))[1], "</text>", sep = ""),
                            paste("<text x=\"", 8.5 * 35.43307 - nchar(names(table(karyotype[,5]))[2]) * 2.1, "\" y=\"", (2.5 + 11.25833 - 0.4) * 35.43307, "\" font-size=\"", karyotype[table(karyotype$species)[1]+1,6], "\" fill=\"#", karyotype[table(karyotype$species)[1]+1,7], "\" transform=\"rotate(-60, 70.86614 487.4999)\" >", names(table(karyotype[,5]))[2], "</text>", sep = ""),
                            paste("<text x=\"", 8.5 * 35.43307 - nchar(names(table(karyotype[,5]))[3]) * 2.1, "\" y=\"", (2.5 + 11.25833 - 0.4) * 35.43307, "\" font-size=\"", karyotype[table(karyotype$species)[1]+table(karyotype$species)[2]+1,6], "\" fill=\"#", karyotype[table(karyotype$species)[1]+table(karyotype$species)[2]+1,7], "\" transform=\"rotate(60, 531.496 487.4999)\" >", names(table(karyotype[,5]))[3], "</text>", sep = ""),
                            stringsAsFactors = FALSE
      )

      for (i in 1 : table(karyotype$species)[1]){
        karyotype[i,15] = paste("<text x=\"", (karyotype[i, 8] + karyotype[i, 8] + karyotype[i, 12])/2 - nchar(karyotype[i, 1]) * 2.1, "\" y=\"", karyotype[i, 9] + 0.4 * 35.43307, "\" font-size=\"4\" fill=\"black\" >", karyotype[i, 1], "</text>", sep = "")
      }
      for (i in (table(karyotype$species)[1] + 1) : (table(karyotype$species)[1] + table(karyotype$species)[2])){
        karyotype[i,15] = paste("<text x=\"", (karyotype[i, 8] + karyotype[i, 8] + karyotype[i, 12])/2 - nchar(karyotype[i, 1]) * 2.1, "\" y=\"", karyotype[i, 9] - 0.1 * 35.43307, "\" font-size=\"4\" fill=\"black\" transform=\"rotate(-60, 70.86614 487.4999)\" >", karyotype[i, 1], "</text>", sep = "")
      }
      for (i in (table(karyotype$species)[1] + table(karyotype$species)[2] + 1) : (table(karyotype$species)[1] + table(karyotype$species)[2] + table(karyotype$species)[3])){
        karyotype[i,15] = paste("<text x=\"", (karyotype[i, 8] + karyotype[i, 8] + karyotype[i, 12])/2 - nchar(karyotype[i, 1]) * 2.1, "\" y=\"", karyotype[i, 9] - 0.1 * 35.43307, "\" font-size=\"4\" fill=\"black\" transform=\"rotate(60, 531.496 487.4999)\" >", karyotype[i, 1], "</text>", sep = "")
      }
      names(karyotype)[15] <- "text"

      gradient1 <- data.frame(paste("<defs>"),
                              paste("<linearGradient id=\"gradient1\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">"),
                              paste("<stop offset=\"0%\" style=\"stop-color:#", karyotype[table(karyotype$species)[1]+1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("<stop offset=\"100%\" style=\"stop-color:#", karyotype[1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("</linearGradient>"),
                              paste("</defs>"),
                              stringsAsFactors = FALSE
      )

      gradient2 <- data.frame(paste("<defs>"),
                              paste("<linearGradient id=\"gradient2\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">"),
                              paste("<stop offset=\"0%\" style=\"stop-color:#", karyotype[1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("<stop offset=\"100%\" style=\"stop-color:#", karyotype[table(karyotype$species)[1]+table(karyotype$species)[2]+1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("</linearGradient>"),
                              paste("</defs>"),
                              stringsAsFactors = FALSE
      )

      gradient3 <- data.frame(paste("<defs>"),
                              paste("<linearGradient id=\"gradient3\" x1=\"0%\" y1=\"0%\" x2=\"100%\" y2=\"0%\">"),
                              paste("<stop offset=\"0%\" style=\"stop-color:#", karyotype[table(karyotype$species)[1]+1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("<stop offset=\"100%\" style=\"stop-color:#", karyotype[table(karyotype$species)[1]+table(karyotype$species)[2]+1,7], ";", sep = ""),
                              paste("stop-opacity:1\"/>"),
                              paste("</linearGradient>"),
                              paste("</defs>"),
                              stringsAsFactors = FALSE
      )

      for (i in 1 : nrow(synteny)) {
        if (synteny[i, 8] == 1 & synteny[i, 7] == "gradient"){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 15], synteny[i, 16], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 9], synteny[i, 10], "L", synteny[i, 11], synteny[i, 12], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 13], synteny[i, 14], "Z\""), " style=\"stroke-width: 1px;stroke:url(#gradient1);fill:url(#gradient1)\"/>", sep = "")
        } else if (synteny[i, 8] == 2 & synteny[i, 7] == "gradient"){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 15], synteny[i, 16], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 9], synteny[i, 10], "L", synteny[i, 11], synteny[i, 12], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 13], synteny[i, 14], "Z\""), " style=\"stroke-width: 1px;stroke:url(#gradient2);fill:url(#gradient2)\"/>", sep = "")
        } else if (synteny[i, 8] == 3 & synteny[i, 7] == "gradient"){
          synteny[i, 21] <- paste(paste("<path d=\"M", synteny[i, 9], synteny[i, 10], "Q", synteny[i, 17], synteny[i, 18], synteny[i, 11], synteny[i, 12], "L", synteny[i, 13], synteny[i, 14], "Q", synteny[i, 19], synteny[i, 20], synteny[i, 15], synteny[i, 16], "Z\""), " style=\"stroke-width: 1px;stroke:url(#gradient3);fill:url(#gradient3)\"/>", sep = "")
        }
      }
    }
  } else {
    mpx<-3.543307 #conversion ratio
    chr_width <- width / (2.6*nrow(karyotype)) * mpx  #width(mm)

    karyotype$x9 <- karyotype$x1 <- karyotype$x8 <- karyotype$x4 <- karyotype$x5 <-
      karyotype$x12 <- apply(data.frame(1:nrow(karyotype)),1,function(x)(20*mpx+(x[1]-1)*2.6*chr_width))

    maxchrlen<-150 #the length of the longest chromosome was set to be 150mm

    karyotype$y1 <- karyotype$y2 <-
      apply(data.frame(karyotype$End),1,
            function(x) ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2)) #the arc radius on each end was set to be chr_width/2

    karyotype$x10 <- karyotype$x2 <- karyotype$x3 <- karyotype$x7 <-
      karyotype$x6 <- karyotype$x11 <- karyotype$x1+chr_width

    if (col_num == 5){
      karyotype$y3 <- karyotype$y8 <-
        apply(data.frame(karyotype$End,karyotype$CE_start),1,
              function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))

      karyotype$y4 <- karyotype$y7 <-
        apply(data.frame(karyotype$End,karyotype$CE_end),1,
              function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))
    }

    karyotype$y5 <- karyotype$y6 <- (25+maxchrlen)*mpx-chr_width/2

    karyotype$y9 <- karyotype$y10 <-
      apply(data.frame(karyotype$End),1,
            function(x)((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx))
    karyotype$y11 <- karyotype$y12 <- (25+maxchrlen)*mpx


    if (col_num == 5){
      karyotype$path = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1,
                             " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2,
                             " L", karyotype$x3, ",", karyotype$y3,
                             " L", karyotype$x4, ",", karyotype$y4,
                             " L", karyotype$x5, ",", karyotype$y5,
                             " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x6, ",", karyotype$y6,
                             " L", karyotype$x7, ",", karyotype$y7,
                             " L", karyotype$x8, ",", karyotype$y8,
                             " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")
    } else {
      karyotype$path = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1,
                             " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2,
                             " L", karyotype$x6, ",", karyotype$y6,
                             " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x5, ",", karyotype$y5,
                             " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")
    }

    karyotype$hat = paste("<path d=\"M", karyotype$x1, ",", karyotype$y1,
                          " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2,
                          " L", karyotype$x10, ",", karyotype$y10,
                          " L", karyotype$x9, ",", karyotype$y9,
                          " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

    karyotype$shoe = paste("<path d=\"M", karyotype$x5, ",", karyotype$y5,
                           " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x6, ",", karyotype$y6,
                           " L", karyotype$x11, ",", karyotype$y11,
                           " L", karyotype$x12, ",", karyotype$y12,
                           " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

    karyotype$bow = paste("<path d=\"M", karyotype$x8, ",", karyotype$y8,
                          " L", karyotype$x7, ",", karyotype$y7,
                          " L", karyotype$x3, ",", karyotype$y3,
                          " L", karyotype$x4, ",", karyotype$y4,
                          " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

    karyotype$text = paste("<text x=\"", (karyotype$x1 + karyotype$x2)/2 -
                             nchar(karyotype$Chr) * 2.2, "\" y=\"",
                           (150 + 25) * mpx + 15,
                           "\" font-size=\"9\" font-family=\"Arial\" fill=\"black\" >",
                           karyotype$Chr, "</text>", sep = "")

    if (!is.null(mydata)) {

      cnum<-10000

      mydata$color <- colorRampPalette(colorset1)(cnum)[round(rescale(mydata$Value,to=c(1,cnum)))]

      mydata<-merge(mydata,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x1=karyotype$x1,x2=karyotype$x2),by="Chr")

      mydata$y1 <- apply(data.frame(mydata$ChrEnd,mydata$Start),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))
      mydata$y2 <- apply(data.frame(mydata$ChrEnd,mydata$End),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))

      mydata$rect = paste("<path d=\"M", mydata$x1, ",", mydata$y1,
                          " L", mydata$x2, ",", mydata$y1,
                          " L", mydata$x2, ",", mydata$y2,
                          " L", mydata$x1, ",", mydata$y2,
                          " Z" ,"\" style=\"fill:", mydata$color, "; stroke:", mydata$color, "; stroke-width:0.25\"/>", sep = "")

      # legend for mydata
      mydata_legend <- data.frame(color=colorRampPalette(colorset1)(cnum))
      mydata_legend$x1<-apply(data.frame(1:cnum),1,function(x)(Lx * mpx + (x[1] - 1) * ((20 * mpx) / cnum)))
      mydata_legend$y1<-Ly * mpx
      mydata_legend$legend<-paste("<rect x=\"", mydata_legend$x1, "\" y=\"", mydata_legend$y1,
                                  "\" width=\"", (20 * mpx) / cnum,
                                  "\" height=\"", 4 * mpx,
                                  "\" style=\"fill:", mydata_legend$color, ";stroke:none\"/>", sep = "")

      legend_text <- c(paste("<text x=\"", min(mydata_legend$x1),
                             "\" y=\"", min(mydata_legend$y1) + 8 * mpx - 3,
                             "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >Low</text>", sep = ""),
                       paste("<text x=\"", max(mydata_legend$x1) - 20,
                             "\" y=\"", min(mydata_legend$y1) + 8 * mpx - 3,
                             "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >High</text>", sep = "")
      )

    }
    if (!is.null(mydata_interval)) {
      if (label_type == "marker"){

        mydata_interval<-merge(mydata_interval,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x2=karyotype$x2),by="Chr")

        mydata_interval$x <- mydata_interval$x2 + chr_width / 2
        mydata_interval$y0 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)((25+maxchrlen*(1-(x[1]-(x[2]+x[3])/2)/max(karyotype$End))) * mpx))
        mydata_interval<-mydata_interval[order(mydata_interval$Chr,mydata_interval$y0),]

        #repel
        repel<- function(mydata,myforce,tag){
          #print(tag)
          #myordata<-sort(mydata)
          #if(length(mydata)==1){return(mydata)}
          if(min(diff(mydata)) >= myforce | tag==3000){return(mydata)}
          tag<-tag+1
          sp<-which.min(diff(mydata))
          ep<-sp+1
          mydata[sp]<-mydata[sp]-myforce
          mydata[ep]<-mydata[ep]+myforce
          mydepo<-sort(mydata)
          return(repel(mydepo,myforce,tag))
        }

        mydata_interval$y <- NA
        for (chr in unique(mydata_interval$Chr)){
          if(nrow(mydata_interval[mydata_interval$Chr==chr,])>1){
            mydata_interval[mydata_interval$Chr == chr,]$y <-repel(mydata_interval[mydata_interval$Chr == chr,]$y0, chr_width/3, 1)
          }else{mydata_interval[mydata_interval$Chr == chr,]$y <- mydata_interval[mydata_interval$Chr == chr,]$y0}
        }

        mydata_interval_triangle<- mydata_interval[mydata_interval$Shape == "triangle",]
        if ("triangle" %in% mydata_interval$Shape){
          mydata_interval_triangle$interval <- paste("<path d=\"M", mydata_interval_triangle$x - chr_width / 4, ",", mydata_interval_triangle$y + chr_width / 4,
                                                     " L", mydata_interval_triangle$x + chr_width / 4, ",", mydata_interval_triangle$y + chr_width / 4,
                                                     " L", mydata_interval_triangle$x, ",", mydata_interval_triangle$y - chr_width / 4,
                                                     " Z" ,"\" style=\"fill:#", mydata_interval_triangle$color, ";stroke:none\"/>", sep = "")
        }

        mydata_interval_box<- mydata_interval[mydata_interval$Shape == "box",]
        if ("box" %in% mydata_interval$Shape){
          mydata_interval_box$interval <- paste("<rect x=\"", mydata_interval_box$x - chr_width / 4,
                                                "\" y=\"", mydata_interval_box$y - chr_width / 4,
                                                "\" width=\"", chr_width / 2,
                                                "\" height=\"", chr_width / 2,
                                                "\" style=\"fill:#", mydata_interval_box$color, "; stroke:none\"/>", sep = "")
        }

        mydata_interval_circle<- mydata_interval[mydata_interval$Shape == "circle",]
        if ("circle" %in% mydata_interval$Shape) {
          mydata_interval_circle$interval <- paste("<circle cx=\"", mydata_interval_circle$x,
                                                   "\" cy=\"", mydata_interval_circle$y,
                                                   "\" r=\"", chr_width / 4,
                                                   "\" style=\"fill:#", mydata_interval_circle$color, "; stroke:none\"/>", sep = "")
        }

        mydata_interval<-rbind(mydata_interval_box,mydata_interval_circle,mydata_interval_triangle)
        mydata_interval$line <- paste("<line x1=\"", mydata_interval$x2,
                                      "\" y1=\"", mydata_interval$y0,
                                      "\" x2=\"", mydata_interval$x,
                                      "\" y2=\"", mydata_interval$y,
                                      "\" style=\"stroke:#", mydata_interval$color, "; stroke-width:0.25\"/>", sep = "")


        # legend for mydata2
        mydata2_legend <- mydata_interval[!duplicated(mydata_interval$Type), 1:6]
        mydata2_legend <- mydata2_legend[order(mydata2_legend$Shape, mydata2_legend$color),]

        #mydata2_legend$x<-mydata_legend$x1[1]
        mydata2_legend$x <- Lx * mpx
        mydata2_legend$y<-apply(data.frame(1:nrow(mydata2_legend)),1,
                                function(x)(Ly * mpx + 4 + (12 + (x[1] - 1) * 4 ) * mpx))
        #function(x)(mydata_legend$y1[1]+ 4 + (12 + (x[1] - 1) * 4 ) * mpx))

        mydata2_legend$x1<-mydata2_legend$x + 4
        mydata2_legend$y1<-mydata2_legend$y- 4 * mpx / 2

        for (i in 1:nrow(mydata2_legend)){
          if (mydata2_legend[i, 3] == "triangle") {
            mydata2_legend[i,11] = paste("<path d=\"M", mydata2_legend[i,9] - 4, ",",
                                         mydata2_legend[i,10] + 4, " L", mydata2_legend[i,9] + 4, ",",
                                         mydata2_legend[i,10] + 4, " L", mydata2_legend[i,9], ",",
                                         mydata2_legend[i,10] - 4, " Z" ,"\" style=\"fill:#",
                                         mydata2_legend[i, 6], "; stroke:none\"/>", sep = "")
          } else if (mydata2_legend[i, 3] == "box") {
            mydata2_legend[i,11] <- paste("<rect x=\"",
                                          mydata2_legend[i,9] - 4, "\" y=\"",
                                          mydata2_legend[i,10] - 4, "\" width=\"", 8,
                                          "\" height=\"", 8, "\" style=\"fill:#",
                                          mydata2_legend[i, 6], ";stroke:none\"/>", sep = "")
          } else if (mydata2_legend[i, 3] == "circle") {
            mydata2_legend[i,11] <- paste("<circle cx=\"",
                                          mydata2_legend[i,9], "\" cy=\"",
                                          mydata2_legend[i,10], "\" r=\"", 4, "\" style=\"fill:#",
                                          mydata2_legend[i, 6], ";stroke:none\"/>", sep = "")
          }
        }
        names(mydata2_legend)[11] <- "shape"

        for (i in 1:nrow(mydata2_legend)){
          mydata2_legend[i,12] <- paste("<text x=\"",
                                        mydata2_legend[i,7] + 15, "\" y=\"",
                                        mydata2_legend[i,8] - (4 * mpx / 2 - 4),
                                        "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >",
                                        mydata2_legend[i,2],"</text>", sep = "")
        }
        names(mydata2_legend)[12] <- "name"
      } else if (label_type == "heatmap") {
        #the positions of chromosome names
        karyotype$text = paste("<text x=\"", (2 * karyotype$x1 + 2.6 * chr_width)/2 - 4 * nchar(karyotype$Chr), "\" y=\"",
                               (150 + 25) * mpx + 15,
                               "\" font-size=\"9\" font-family=\"Arial\" fill=\"black\" >",
                               karyotype$Chr, "</text>", sep = "")

        #draw new idiograms
        if (col_num == 5){
          karyotype$path2 = paste("<path d=\"M", karyotype$x1 + 1.2 * chr_width, ",", karyotype$y1,
                                  " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2 + 1.2 * chr_width, ",", karyotype$y2,
                                  " L", karyotype$x3 + 1.2 * chr_width, ",", karyotype$y3,
                                  " L", karyotype$x4 + 1.2 * chr_width, ",", karyotype$y4,
                                  " L", karyotype$x5 + 1.2 * chr_width, ",", karyotype$y5,
                                  " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x6 + 1.2 * chr_width, ",", karyotype$y6,
                                  " L", karyotype$x7 + 1.2 * chr_width, ",", karyotype$y7,
                                  " L", karyotype$x8 + 1.2 * chr_width, ",", karyotype$y8,
                                  " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")
        } else {
          karyotype$path2 = paste("<path d=\"M", karyotype$x1 + 1.2 * chr_width, ",", karyotype$y1,
                                  " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2 + 1.2 * chr_width, ",", karyotype$y2,
                                  " L", karyotype$x6 + 1.2 * chr_width, ",", karyotype$y6,
                                  " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x5 + 1.2 * chr_width, ",", karyotype$y5,
                                  " Z" ,"\" style=\"fill:none; stroke:grey; stroke-width:1\"/>", sep = "")
        }

        karyotype$hat2 = paste("<path d=\"M", karyotype$x1 + 1.2 * chr_width, ",", karyotype$y1,
                               " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2 + 1.2 * chr_width, ",", karyotype$y2,
                               " L", karyotype$x10 + 1.2 * chr_width, ",", karyotype$y10,
                               " L", karyotype$x9 + 1.2 * chr_width, ",", karyotype$y9,
                               " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

        karyotype$shoe2 = paste("<path d=\"M", karyotype$x5 + 1.2 * chr_width, ",", karyotype$y5,
                                " A", chr_width/2, ",", chr_width/2, " 0 0,0 ", karyotype$x6 + 1.2 * chr_width, ",", karyotype$y6,
                                " L", karyotype$x11 + 1.2 * chr_width, ",", karyotype$y11,
                                " L", karyotype$x12 + 1.2 * chr_width, ",", karyotype$y12,
                                " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")

        karyotype$bow2 = paste("<path d=\"M", karyotype$x8 + 1.2 * chr_width, ",", karyotype$y8,
                               " L", karyotype$x7 + 1.2 * chr_width, ",", karyotype$y7,
                               " L", karyotype$x3 + 1.2 * chr_width, ",", karyotype$y3,
                               " L", karyotype$x4 + 1.2 * chr_width, ",", karyotype$y4,
                               " Z" ,"\" style=\"fill:white; stroke:white; stroke-width:0.75\"/>", sep = "")


        cnum<-10000

        mydata_interval$color <- colorRampPalette(colorset2)(cnum)[round(rescale(mydata_interval$Value,to=c(1,cnum)))]

        mydata_interval<-merge(mydata_interval,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x1=karyotype$x1,x2=karyotype$x2),by="Chr")

        mydata_interval$y1 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))
        mydata_interval$y2 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$End),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))

        mydata_interval$rect = paste("<path d=\"M", mydata_interval$x1 + 1.2 * chr_width, ",", mydata_interval$y1,
                                     " L", mydata_interval$x2 + 1.2 * chr_width, ",", mydata_interval$y1,
                                     " L", mydata_interval$x2 + 1.2 * chr_width, ",", mydata_interval$y2,
                                     " L", mydata_interval$x1 + 1.2 * chr_width, ",", mydata_interval$y2,
                                     " Z" ,"\" style=\"fill:", mydata_interval$color, "; stroke:", mydata_interval$color, "; stroke-width:0.25\"/>", sep = "")

        # legend for mydata_interval
        mydata2_legend <- data.frame(color=colorRampPalette(colorset2)(cnum))
        mydata2_legend$x1<-apply(data.frame(1:cnum),1,function(x)(Lx * mpx + (x[1] - 1) * ((20 * mpx) / cnum)))
        mydata2_legend$y1<-Ly * mpx + 4 + 12 * mpx
        mydata2_legend$legend<-paste("<rect x=\"", mydata2_legend$x1, "\" y=\"", mydata2_legend$y1,
                                     "\" width=\"", (20 * mpx) / cnum,
                                     "\" height=\"", 4 * mpx,
                                     "\" style=\"fill:", mydata2_legend$color, ";stroke:none\"/>", sep = "")

        legend_text2 <- c(paste("<text x=\"", min(mydata2_legend$x1),
                                "\" y=\"", min(mydata2_legend$y1) + 8 * mpx - 3,
                                "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >Low</text>", sep = ""),
                          paste("<text x=\"", max(mydata2_legend$x1) - 20,
                                "\" y=\"", min(mydata2_legend$y1) + 8 * mpx - 3,
                                "\" font-size=\"12\" font-family=\"Arial\" fill=\"black\" >High</text>", sep = "")
        )

      } else if (label_type == "line"){
        if (ncol(mydata_interval) == 5){
          mydata_interval <- merge(mydata_interval, data.frame(Chr = karyotype$Chr, ChrEnd = karyotype$End, x = karyotype$x1), by="Chr")
          mydata_interval$x <- mydata_interval$x + 1.2 * chr_width + mydata_interval$Value * chr_width / max(mydata_interval$Value)
          mydata_interval$y <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)(((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx + (25+maxchrlen*(1-(x[1]-x[3])/max(karyotype$End))) * mpx)/2))
          mydata_interval$point <- paste(mydata_interval$x, mydata_interval$y, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 29] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point[j]
              karyotype[i, 29] <- paste(karyotype[i, 29], tmp, sep = " ")
            }
            karyotype[i, 29] <- paste("<polyline points=\"",
                                      karyotype[i, 29],
                                      "\" style=\"fill:none; stroke:#", mydata_interval[1,5], "; stroke-width:1\"/>", sep = "")
          }
          colnames(karyotype)[29] <- "line"
        } else if (ncol(mydata_interval) == 7){
          mydata_interval <- merge(mydata_interval, data.frame(Chr = karyotype$Chr, ChrEnd = karyotype$End, x = karyotype$x1), by="Chr")
          mydata_interval$x1 <- mydata_interval$x + 1.2 * chr_width + mydata_interval$Value_1 * chr_width / max(mydata_interval$Value_1)
          mydata_interval$y1 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)(((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx + (25+maxchrlen*(1-(x[1]-x[3])/max(karyotype$End))) * mpx)/2))
          mydata_interval$point1 <- paste(mydata_interval$x1, mydata_interval$y1, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 29] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point1[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point1[j]
              karyotype[i, 29] <- paste(karyotype[i, 29], tmp, sep = " ")
            }
            karyotype[i, 29] <- paste("<polyline points=\"",
                                      karyotype[i, 29],
                                      "\" style=\"fill:none; stroke:#", mydata_interval[1,5], "; stroke-width:1\"/>", sep = "")
          }
          colnames(karyotype)[29] <- "line1"

          mydata_interval$x2 <- mydata_interval$x + 1.2 * chr_width + mydata_interval$Value_2 * chr_width / max(mydata_interval$Value_2)
          mydata_interval$y2 <- mydata_interval$y1
          mydata_interval$point2 <- paste(mydata_interval$x2, mydata_interval$y2, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 30] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point2[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point2[j]
              karyotype[i, 30] <- paste(karyotype[i, 30], tmp, sep = " ")
            }
            karyotype[i, 30] <- paste("<polyline points=\"",
                                      karyotype[i, 30],
                                      "\" style=\"fill:none; stroke:#", mydata_interval[1,7], "; stroke-width:1\"/>", sep = "")
          }
          colnames(karyotype)[30] <- "line2"
        }
      } else if (label_type == "polygon"){
        if (ncol(mydata_interval) == 5){
          mydata_interval <- merge(mydata_interval, data.frame(Chr = karyotype$Chr, ChrEnd = karyotype$End, x = karyotype$x1), by="Chr")
          mydata_interval$x <- mydata_interval$x + 1.2 * chr_width + mydata_interval$Value * chr_width / max(mydata_interval$Value)
          mydata_interval$y <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)(((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx + (25+maxchrlen*(1-(x[1]-x[3])/max(karyotype$End))) * mpx)/2))
          mydata_interval$point <- paste(mydata_interval$x, mydata_interval$y, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 29] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point[j]
              karyotype[i, 29] <- paste(karyotype[i, 29], tmp, sep = " ")
            }
            karyotype[i, 29] <- paste("<polygon points=\"",
                                      paste(karyotype$x1[i] + 1.2 * chr_width, ",", karyotype$y1[i] - chr_width/2, " ", sep = ""),
                                      karyotype[i, 29],
                                      paste(" ", karyotype$x1[i] + 1.2 * chr_width, ",", (25 + maxchrlen) * mpx, sep = ""),
                                      "\" style=\"fill:#", mydata_interval[1,5], "; stroke:#", mydata_interval[1,5], "; stroke-width:0.25\"/>", sep = "")
          }
          colnames(karyotype)[29] <- "polygon"
        } else if (ncol(mydata_interval) == 7){
          mydata_interval <- merge(mydata_interval, data.frame(Chr = karyotype$Chr, ChrEnd = karyotype$End, x = karyotype$x1), by="Chr")
          mydata_interval$x1 <- mydata_interval$x + 1.2 * chr_width + mydata_interval$Value_1 * chr_width / max(mydata_interval$Value_1)
          mydata_interval$y1 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)(((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx + (25+maxchrlen*(1-(x[1]-x[3])/max(karyotype$End))) * mpx)/2))
          mydata_interval$point1 <- paste(mydata_interval$x1, mydata_interval$y1, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 29] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point1[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point1[j]
              karyotype[i, 29] <- paste(karyotype[i, 29], tmp, sep = " ")
            }
            karyotype[i, 29] <- paste("<polygon points=\"",
                                      paste(karyotype$x1[i] + 1.2 * chr_width, ",", karyotype$y1[i] - chr_width/2, " ", sep = ""),
                                      karyotype[i, 29],
                                      paste(" ", karyotype$x1[i] + 1.2 * chr_width, ",", (25 + maxchrlen) * mpx, sep = ""),
                                      "\" style=\"fill:#", mydata_interval[1,5], "; stroke:#", mydata_interval[1,5], "; stroke-width:0.25\"/>", sep = "")
          }
          colnames(karyotype)[29] <- "polygon1"

          mydata_interval$x2 <- mydata_interval$x + 1.4 * chr_width + mydata_interval$Value_2 * chr_width / max(mydata_interval$Value_2)
          mydata_interval$y2 <- mydata_interval$y1
          mydata_interval$point2 <- paste(mydata_interval$x2, mydata_interval$y2, sep = ",")
          for (i in 1:nrow(karyotype)){
            karyotype[i, 30] <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point2[1]
            for (j in 2:nrow(subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1]))){
              tmp <- subset(mydata_interval, mydata_interval$Chr == karyotype[i, 1])$point2[j]
              karyotype[i, 30] <- paste(karyotype[i, 30], tmp, sep = " ")
            }
            karyotype[i, 30] <- paste("<polygon points=\"",
                                      paste(karyotype$x1[i] + 1.4 * chr_width, ",", karyotype$y1[i] - chr_width/2, " ", sep = ""),
                                      karyotype[i, 30],
                                      paste(" ", karyotype$x1[i] + 1.4 * chr_width, ",", (25 + maxchrlen) * mpx, sep = ""),
                                      "\" style=\"fill:#", mydata_interval[1,7], "; stroke:#", mydata_interval[1,7], "; stroke-width:0.25\"/>", sep = "")
          }
          colnames(karyotype)[30] <- "polygon2"
        }
      }
    }
  }

  #write the head
  first_line <- c("<?xml version=\"1.0\" standalone=\"no\"?>",
                  "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                  "",
                  paste("<svg id=\"svg\" width=\"744.0945\" height=\"1052.362\">", "\t")
  )


  cat(first_line, file = output)

  if (!is.null(synteny)){
    if (length(table(karyotype$species)) == 2){
      cat(species[1, 1], file = output, append = TRUE)
      cat(species[1, 2], file = output, append = TRUE)
      cat(karyotype$text, file = output, append = TRUE)
      cat(synteny$path, file = output, append = TRUE)

      cat(karyotype$rect_outer, file = output, append = TRUE)
      cat(karyotype$rect_inner, file = output, append = TRUE)
    } else if (length(table(karyotype$species)) == 3) {
      cat(gradient1[1, 1], file = output, append = TRUE)
      cat(gradient1[1, 2], file = output, append = TRUE)
      cat(gradient1[1, 3], file = output, append = TRUE)
      cat(gradient1[1, 4], file = output, append = TRUE)
      cat(gradient1[1, 5], file = output, append = TRUE)
      cat(gradient1[1, 6], file = output, append = TRUE)
      cat(gradient1[1, 7], file = output, append = TRUE)
      cat(gradient1[1, 8], file = output, append = TRUE)

      cat(gradient2[1, 1], file = output, append = TRUE)
      cat(gradient2[1, 2], file = output, append = TRUE)
      cat(gradient2[1, 3], file = output, append = TRUE)
      cat(gradient2[1, 4], file = output, append = TRUE)
      cat(gradient2[1, 5], file = output, append = TRUE)
      cat(gradient2[1, 6], file = output, append = TRUE)
      cat(gradient2[1, 7], file = output, append = TRUE)
      cat(gradient2[1, 8], file = output, append = TRUE)

      cat(gradient3[1, 1], file = output, append = TRUE)
      cat(gradient3[1, 2], file = output, append = TRUE)
      cat(gradient3[1, 3], file = output, append = TRUE)
      cat(gradient3[1, 4], file = output, append = TRUE)
      cat(gradient3[1, 5], file = output, append = TRUE)
      cat(gradient3[1, 6], file = output, append = TRUE)
      cat(gradient3[1, 7], file = output, append = TRUE)
      cat(gradient3[1, 8], file = output, append = TRUE)

      cat(species[1, 1], file = output, append = TRUE)
      cat(species[1, 2], file = output, append = TRUE)
      cat(species[1, 3], file = output, append = TRUE)
      cat(karyotype$text, file = output, append = TRUE)
      cat(synteny$path, file = output, append = TRUE)

      cat(karyotype$rect_outer, file = output, append = TRUE)
      cat(karyotype$rect_inner, file = output, append = TRUE)
    }
  } else {
    cat(mydata$rect, file = output, append = TRUE)

    #write the idiograms
    cat(karyotype$hat, file = output, append = TRUE)
    cat(karyotype$shoe, file = output, append = TRUE)
    cat(karyotype$bow, file = output, append = TRUE)
    cat(karyotype$path, file = output, append = TRUE)
    cat(karyotype$text, file = output, append = TRUE)

    #legend
    ## write.table(mydata_legend$legend, output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    ## write.table(legend_text[1,1], output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
    ## write.table(legend_text[1,2], output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

    if (!is.null(mydata)) {

      cat(legend_text, file = output, append = TRUE)
      cat(mydata_legend$legend, file = output, append = TRUE)

    }

    cat(mydata_interval$interval, file = output, append = TRUE)
    cat(mydata_interval$rect, file = output, append = TRUE)

    if (!is.null(mydata_interval)) {
      if (label_type == "marker"){

        cat(mydata_interval$line, file = output, append = TRUE)

        cat(mydata2_legend$shape, file = output, append = TRUE)
        cat(mydata2_legend$name, file = output, append = TRUE)

      } else if (label_type == "heatmap") {
        cat(karyotype$hat2, file = output, append = TRUE)
        cat(karyotype$shoe2, file = output, append = TRUE)
        cat(karyotype$bow2, file = output, append = TRUE)
        cat(karyotype$path2, file = output, append = TRUE)
        cat(legend_text2, file = output, append = TRUE)
        cat(mydata2_legend$legend, file = output, append = TRUE)
      } else if (label_type == "line") {
        if (ncol(mydata_interval) == 9){
          cat(karyotype$line, file = output, append = TRUE)
        } else if (ncol(mydata_interval) == 15){
          cat(karyotype$line1, file = output, append = TRUE)
          cat(karyotype$line2, file = output, append = TRUE)
        }
      } else if (label_type == "polygon") {
        if (ncol(mydata_interval) == 9){
          cat(karyotype$polygon, file = output, append = TRUE)
        } else if (ncol(mydata_interval) == 15){
          cat(karyotype$polygon2, file = output, append = TRUE)
          cat(karyotype$polygon1, file = output, append = TRUE)
        }
      }
    }
  }

  ##write the tail
  cat("</svg>", file = output, append = TRUE)

}
