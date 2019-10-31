##' extract some specific feature information from a gff file
##'
##'
##' @title GFFex
##' @param input gff file
##' @param karyotype karyotype file
##' @param feature feature name
##' @param window window size
##' @return dataframe
##' @importFrom tidyr separate
##' @export
##' @rdname GFFex
##' @author Zhaodong Hao, Dekang Lv, Ying Ge, Jisen Shi, Dolf Weijers, Guangchuang Yu, Jinhui Chen

GFFex = function(input, karyotype, feature = "gene", window = 1000000){
  gff <- read.table(input,
                    stringsAsFactors = F,
                    header = F,
                    comment.char = "#",
                    sep = '\t',
                    quote = ""
                    )
  karyotype <- read.table(karyotype,
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F
                          )
  gff <- subset(gff, V1 %in% karyotype$Chr & V3 == feature)

  list_chr <- vector("list", length(names(table(gff$V1))))
  names(list_chr) <- names(table(gff$V1))
  for (i in 1:(length(list_chr))){
    list_chr[[i]] <- as.data.frame(table(cut(subset(gff, V1 == names(list_chr[i]))$V4,
                                           breaks = c(seq(0, subset(karyotype, Chr == names(list_chr[i]))[1,3], window),
                                                      subset(karyotype, Chr == names(list_chr[i]))[1,3]))))
    list_chr[[i]] <- tidyr::separate(list_chr[[i]], Var1, into = c("Start","End"), sep = ",")
    list_chr[[i]]$Start <- gsub('\\(', '', list_chr[[i]]$Start)
    list_chr[[i]]$End <- gsub('\\]', '', list_chr[[i]]$End)
    list_chr[[i]]$Start <- as.numeric(list_chr[[i]]$Start)
    list_chr[[i]]$End <- as.numeric(list_chr[[i]]$End)
    list_chr[[i]]$Start <- list_chr[[i]]$Start + 1
    list_chr[[i]][nrow(list_chr[[i]]),2] <- subset(karyotype, Chr == names(list_chr[i]))[1,3]
    list_chr[[i]]$Chr <- names(list_chr[i])
    list_chr[[i]] <- cbind(list_chr[[i]][,4], list_chr[[i]][,1:3])
    colnames(list_chr[[i]]) <- c("Chr", "Start", "End", "Value")
    list_chr[[i]]$Chr <- as.character(list_chr[[i]]$Chr)
  }

  l <- data.frame()
  for (i in 1:(length(list_chr))){
    df.now <- list_chr[[i]]
    l <- rbind(l, df.now)
  }

  data.frame(l)
}

