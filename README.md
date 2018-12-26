# RIdeogram

RIdeogram is a R package to draw SVG graphics to visualize and map genome-wide data in idiograms (Scalable Vector Graphics http://tutorials.jenkov.com/svg/index.html). The scripts were first developed to draw SVG graphics and Guangchuang Yu helped to integrate these scripts into a R package (https://github.com/GuangchuangYu/ideogram).

# Installation

To install package from Github, you need to install the "devtools" package first<br>
```
install.packages("devtools")
library(devtools)
```
Then, you can install the package "RIdeogram"<br>
```
devtools::install_github('TickingClock1992/RIdeogram')
```

# Citation


Zhaodong Hao, Dekang Lv, Guangchuang Yu, Ying Ge, Jisen Shi and Jinhui
Chen (2018). RIdeogram: Drawing SVG graphics to visualize and map
genome-wide data in idiograms. R package version 0.1.0.

# Usage and Examples

This is a simple package with only two functions `ideogram` and `convertSVG`.

First, you need to load the package after you installed it.

```
require(RIdeogram)
```
Then, you need to load the data from the RIdeogram package.
```
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")
```
You can use the function "head()" to see the data format.
```
head(human_karyotype)
```
Specifically, the 'karyotype' file contains the karyotype information and has five columns (or three, see below). The first column is Chromosome ID, the second and thrid columns are start and end positions of corresponding chromosomes and the fourth and fifth columns are start and end positions of corresponding centromeres.


```
head(gene_density)
```
The 'mydata' file contains the heatmap information and has four columns. The first column is Chromosome ID, the second and thrid columns are start and end positions of windows in corresponding chromosomes and the fourth column is a characteristic value in corresponding windows, such as gene number.


```
head(Random_RNAs_500)
```
The 'mydata_interval' file contains the label information and has six columns. The first column is the label type, the second column is the shape of label with three available options of box, triangle and circle, the third column is Chromosome ID, the fourth and fifth columns are the start and end positions of corresponding labels in the chromosomes and the sixth column is the color of the label.

Or, you can also load your own data by using the function "read.table", such as
```
human_karyotype <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
gene_density <- read.table("data_1.txt", sep = "\t", header = T, stringsAsFactors = F)
Random_RNAs_500 <- read.table("data_2.txt", sep = "\t", header = T, stringsAsFactors = F)
```
The "karyotype.txt" file contains karyotype information; the "data_1.txt" file contains heatmap data; the "data_2.txt" contains track label data.

These three files are all you need, now you can visualize these information using the 'ideogram' function.

Basic usage
```
ideogram(karyotype, overlaid = NULL, label = NULL, colorset1, colorset2, width, Lx, Ly, output = "chromosome.svg")
convertSVG(svg, device, width, height, dpi)
```

Now, let's begin.

First, we draw a idiogram with no mapping data.
```
ideogram(karyotype = human_karyotype)
convertSVG("chromosome.svg", device = "png")
```
Then, you will find a SVG file and a PNG file in your Working Directory.

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example1.png)

Next, we can map genome-wide data on the chromosome idiogram. In this case, we visulize the gene density across the human genome.
```
ideogram(karyotype = human_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example2.png)

Alternatively, we can map some genome-wide data with track labels next to the chromosome idiograms.
```
ideogram(karyotype = human_karyotype, label = Random_RNAs_500)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example3.png)

We can also map the overlaid heatmap and track labels on the chromosome idiograms at the same time.
```
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example4.png)

If you want to change the color of heatmap, you can modify the argument 'colorset' (default set is colorset = c("#4575b4", "#ffffbf", "#d73027")). You can use either color names as listed by `colors()` or hexadecimal strings of the form "#rrggbb" or "#rrggbbaa".

```
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, colorset = c("#fc8d59", "#ffffbf", "#91bfdb"))
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example5.png)

If you don not know the centromere information in your species, you don not need to modify the script. In this case, the 'karyotype' file has only three columns.<br>
To simulate this case, we deleted the last two columns of the 'human_karyotype' file.
```
human_karyotype <- human_karyotype[,1:3]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500)
convertSVG("chromosome.svg", device = "png")
```
![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example6.png)

If there are only ten chromosomes in your species, maybe you need to motify the argument 'width' (default value is "170").

To simulate this case, we only keep the first ten columns of the 'human_karyotype' file.

Before
```
human_karyotype <- human_karyotype[1:10,]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example7.png)

After
```
human_karyotype <- human_karyotype[1:10,]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, width = 100)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example8.png)

If you want to move the Legend, then you need to modify the arguments 'Lx' and 'Ly'(default values are "160" and "35", separately).<br>
'Lx' means the distance between upper-left point of the Legend and the leaf margin; 'Ly' means the distance between upper-left point of the Legend and the upper margin.

```
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, width = 100, Lx = 80, Ly = 25)
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example9.png)

If you have two sets of heatmap data, such as gene density and LTR density, you can use the following scripts to map and visualize these data in idiograms.
```
data(human_karyotype, package="RIdeogram") #reload the karyotype data
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = LTR_density, colorset1 = c("#f7f7f7", "#e34a33"), colorset2 = c("#f7f7f7", "#2c7fb8")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
convertSVG("chromosome.svg", device = "png")
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/image/example10.png)

In addition, you can use the argument "device" (default value is "png")to set the format of output file, such as, "tiff", "pdf", "jpg", etc. And, you can use the argument "dpi" (default value is "300") to set the resolution of the output image file.
```
convertSVG("chromosome.svg", device = "tiff", dpi = 600)
```
# THANKS
Welcome to any suggestions and discussions.
