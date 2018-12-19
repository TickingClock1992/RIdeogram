##' ideogram with overlaid heatmap annotation and optional track label
##'
##'
##' @title ideogram
##' @param karyotype karyotype data
##' @param overlaid overlaid data
##' @param label track label data
##' @param colorset overlaid heatmap color
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
##' ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500)
##' svg2jpg("chromosome.svg", "chromosome.jpg")
##'
##' # Then, you will find a SVG file and a JPEG file in your Working Directory.
##'
##' @author Zhaodong Hao, Dekang Lv, Guangchuang Yu, Ying Ge, Jisen Shi, Jinhui Chen
##'
ideogram <- function(karyotype, overlaid = NULL, label = NULL, colorset = c("#4575b4", "#ffffbf", "#d73027"), width = 170, Lx = 160, Ly = 35, output = "chromosome.svg") {
  karyotype <- karyotype
  mydata <- overlaid
  mydata_interval <- label

  col_num <- ncol(karyotype)

  mpx<-3.543307 #换算比
  chr_width <- width / (2*nrow(karyotype)) * mpx  #width(mm)可以根据染色体条数适当调整

  karyotype$x9 <- karyotype$x1 <- karyotype$x8 <- karyotype$x4 <- karyotype$x5 <-
    karyotype$x12 <- apply(data.frame(1:nrow(karyotype)),1,function(x)(20*mpx+(x[1]-1)*2*chr_width))

  maxchrlen<-150 #最长染色体长度设置为150mm

  karyotype$y1 <- karyotype$y2 <-
    apply(data.frame(karyotype$End),1,
          function(x) ((25+maxchrlen*(1-x[1]/max(karyotype$End))) * mpx + chr_width/2)) #染色体顶部圆弧高度为chr_width/2

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
                           " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2, #圆弧
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
                        " A", chr_width/2, ",", chr_width/2, " 0 1,1 ", karyotype$x2, ",", karyotype$y2, #圆弧
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

    cnum<-10000 #颜色板色彩的数量设置

    mydata$color <- colorRampPalette(colorset)(cnum)[round(rescale(mydata$Value,to=c(1,cnum)))] #根据value配置颜色，使用白色则要求mydata中基因组均有value

    mydata<-merge(mydata,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x1=karyotype$x1,x2=karyotype$x2),by="Chr")

    mydata$y1 <- apply(data.frame(mydata$ChrEnd,mydata$Start),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))
    mydata$y2 <- apply(data.frame(mydata$ChrEnd,mydata$End),1,function(x)((25+maxchrlen*(1-(x[1]-x[2])/max(karyotype$End))) * mpx))

    mydata$rect = paste("<path d=\"M", mydata$x1, ",", mydata$y1,
                        " L", mydata$x2, ",", mydata$y1,
                        " L", mydata$x2, ",", mydata$y2,
                        " L", mydata$x1, ",", mydata$y2,
                        " Z" ,"\" style=\"fill:", mydata$color, "; stroke:", mydata$color, "; stroke-width:0.25\"/>", sep = "")

    # legend for mydata
    mydata_legend <- data.frame(color=colorRampPalette(colorset)(cnum)) #颜色与上面data1的设置保持一致
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

    mydata_interval<-merge(mydata_interval,data.frame(Chr=karyotype$Chr,ChrEnd=karyotype$End,x2=karyotype$x2),by="Chr")

    mydata_interval$x <- mydata_interval$x2 + chr_width / 3
    mydata_interval$y0 <- apply(data.frame(mydata_interval$ChrEnd,mydata_interval$Start,mydata_interval$End),1,function(x)((25+maxchrlen*(1-(x[1]-(x[2]+x[3])/2)/max(karyotype$End))) * mpx))
    mydata_interval<-mydata_interval[order(mydata_interval$Chr,mydata_interval$y0),]

    #repel函数，用于移动染色体旁边的标记，避免太多时互相重叠
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
    mydata_interval_triangle$interval <- paste("<path d=\"M", mydata_interval_triangle$x - chr_width / 6, ",", mydata_interval_triangle$y + chr_width / 6,
                                               " L", mydata_interval_triangle$x + chr_width / 6, ",", mydata_interval_triangle$y + chr_width / 6,
                                               " L", mydata_interval_triangle$x, ",", mydata_interval_triangle$y - chr_width / 6,
                                               " Z" ,"\" style=\"fill:#", mydata_interval_triangle$color, ";stroke:none\"/>", sep = "")

    mydata_interval_box<- mydata_interval[mydata_interval$Shape == "box",]
    mydata_interval_box$interval <- paste("<rect x=\"", mydata_interval_box$x - chr_width / 6,
                                          "\" y=\"", mydata_interval_box$y - chr_width / 6,
                                          "\" width=\"", chr_width / 3,
                                          "\" height=\"", chr_width / 3,
                                          "\" style=\"fill:#", mydata_interval_box$color, "; stroke:none\"/>", sep = "")

    mydata_interval_circle<- mydata_interval[mydata_interval$Shape == "circle",]
    mydata_interval_circle$interval <- paste("<circle cx=\"", mydata_interval_circle$x,
                                             "\" cy=\"", mydata_interval_circle$y,
                                             "\" r=\"", chr_width / 6,
                                             "\" style=\"fill:#", mydata_interval_circle$color, "; stroke:none\"/>", sep = "")

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
  }


  #写入文件头
  first_line <- c("<?xml version=\"1.0\" standalone=\"no\"?>",
                  "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"",
                  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">",
                  "",
                  paste("<svg id=\"svg\" width=\"744.0945\" height=\"1052.362\">", "\t")
  )


  cat(first_line, file = output)

  cat(mydata$rect, file = output, append = TRUE)

  #染色体轮廓
  cat(karyotype$hat, file = output, append = TRUE)
  cat(karyotype$shoe, file = output, append = TRUE)
  cat(karyotype$bow, file = output, append = TRUE)
  cat(karyotype$path, file = output, append = TRUE)
  cat(karyotype$text, file = output, append = TRUE)

  #色彩legend
  ## write.table(mydata_legend$legend, output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  ## write.table(legend_text[1,1], output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  ## write.table(legend_text[1,2], output, col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

  if (!is.null(mydata)) {

    cat(legend_text, file = output, append = TRUE)
    cat(mydata_legend$legend, file = output, append = TRUE)

  }

  cat(mydata_interval$interval, file = output, append = TRUE)

  if (!is.null(mydata_interval)) {
    ## 连线，不想要染色体和标记之间连线的，就不运行这行
    cat(mydata_interval$line, file = output, append = TRUE)

    ## 形状和形状的legend，不想要染色体旁边标记的，不运行这部分
    cat(mydata2_legend$shape, file = output, append = TRUE)
    cat(mydata2_legend$name, file = output, append = TRUE)

  }


  ## 必备的语句
  ##文件尾
  cat("</svg>", file = output, append = TRUE)

}
