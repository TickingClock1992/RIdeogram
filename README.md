# RIdeogram

RIdeogram is a R package to draw the chromosome ideograms based on Scalable Vector Graphics (SVG; http://tutorials.jenkov.com/svg/index.html).

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

None, for now.

# Usage and Examples

This is a simple package wiht only one function 'ideogram'.<br>

First, you need to load the package after you installed it.
```
require(RIdeogram)
```
Then, you need to load the data. <br>
(Or, you can load your own data by using the function "read.table")
```
data(karyotype, package="RIdeogram")
data(mydata, package="RIdeogram")
data(mydata_interval, package="RIdeogram")
```
You can use the function "head()" to see the data format.
```
head(karyotype)
head(mydata)
head(mydata_interval)
```
Specifically, the 'karyotype' file contains the karyotype information and has five columns (or three, see below). The first column is Chromosome ID, the second and thrid columns are start and end positions of corresponding chromosomes and the fourth and fifth columns are start and end positions of corresponding centromeres.<br>

The 'mydata' file contains the heatmap information and has four columns. The first column is Chromosome ID, the second and thrid columns are start and end positions of windows in corresponding chromosomes and the fourth column is a characteristic value in corresponding windows, such as gene number.<br>

The 'mydata_interval' file contains the label information and has six columns. The first column is the label type, the second column is the shape of label with three available options of box, triangle and circle, the third column is Chromosome ID, the fourth and fifth columns are the start and end positions of corresponding labels in the chromosomes and the sixth column is the color of the label.<br>

These three files are all you need, now you can visualize these information using the 'ideogram' function.<br>
Basic usage
```
ideogram(karyotype, overlaid, label = NULL, colorset, width, Lx, Ly, output = "chromosome.svg")
```

Now, let's begin.
```
ideogram(karyotype, mydata, mydata_interval, c("navy", "white", "firebrick3"), 170, 160, 35, svgfile)
svg2pdf(svgfile, pdffile)
```
Then, you will find a PDF file in your Working Directory.

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example1.jpg)

If you want to change the color of heatmap, you can modify the argument 'colorset', i.e., "c("navy", "white", "firebrick3")". You can use either color names as listed by `colors()` or hexadecimal strings of the form "#rrggbb" or "#rrggbbaa".<br>
```
ideogram(karyotype, mydata, mydata_interval, c("#4E7DB8", "#FEF8B5", "#D73027"), 170, 160, 35, svgfile)
svg2pdf(svgfile, pdffile)
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example2.jpg)

If you don not know the centromere information in your species, you don not need to modify the script. In this case, the 'karyotype' file has only three columns.<br>
To simulate this case, we deleted the last two columns of the 'karyotype' file.
```
karyotype <- karyotype[,1:3]
ideogram(karyotype, mydata, mydata_interval, c("navy", "white", "firebrick3"), 170, 160, 35, svgfile)
svg2pdf(svgfile, pdffile)
```
![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example3.jpg)

If there are only ten chromosomes in your sepcies, maybe you need to motify the argument 'width', i.e., "170".<br>
To simulate this case, we only keep the first ten columns of the 'karyotype' file.<br>

Before
```
karyotype <- karyotype[1:10,]
ideogram(karyotype, mydata, mydata_interval, c("navy", "white", "firebrick3"), 170, 160, 35, svgfile)
svg2pdf(svgfile, pdffile)
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example4.jpg)

After
```
karyotype <- karyotype[1:10,]
ideogram(karyotype, mydata, mydata_interval, c("navy", "white", "firebrick3"), 100, 160, 35, svgfile)
svg2pdf(svgfile, pdffile)
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example5.jpg)

If you want to move the Legend, then you need to modify the argument 'Lx' and 'Ly', i.e., "160" and "35".<br>
'Lx' means the distance between upper-left point of the Legend and the leaf margin; 'Ly' means the distance between upper-left point of the Legend and the upper margin.

```
karyotype <- karyotype[1:10,]
ideogram(karyotype, mydata, mydata_interval, c("navy", "white", "firebrick3"), 100, 80, 25, svgfile)
svg2pdf(svgfile, pdffile)
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example6.jpg)

If you don't want to plot the label, then you just need to modify the argument 'label', i.e., "mydata_interval".
```
karyotype <- karyotype[1:10,]
ideogram(karyotype, mydata, NULL, c("navy", "white", "firebrick3"), 100, 80, 25, svgfile)
svg2pdf(svgfile, pdffile)
```

![image](https://github.com/TickingClock1992/RIdeogram/blob/master/images/example7.jpg)

# THANKS
Welcome to any suggestions and discussions.
