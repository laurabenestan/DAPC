# Discriminant Analysis in Principal Components (DPCA)

This tutorial provides recommendations for applying the Discriminant Analysis of Principal Components (DAPC) to a wide genomic dataset such as RAD-seq or DarT. 

DPCA acts like a Principal Component Analysis (PCA) and will give more importance to intergroup differences rather than the intragroup differences. Groups can be specify at the beginning of the analysis in a blind way (hereafter named DPCA with no prior) or based on sampling information (DPCA with prior). 

DPCA analysis is particularly relevant in a case of marine species, where low genetic differentiation is often found and hamper the possibility to reveal substle but significant genetic differences.

This statistical approach was developped by Jombart and colleagues through the `adegenet` package. See the [excellent tutorial](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf) done by Thibault Jombart and Cattlin Collins.

## 0 - Prepare your R environment

- [R Version 3.6.0](https://cran.r-project.org/)
	* R packages: adegenet, ade4, ggplot2, reshape, stringr, radiator
  
## 1 - Download genomic and spatial datasets

Transform the vcf to a genepop file using `radiator` package.
```{r}
snps_genlight <- genomic_converter(data = "4991_515ind.vcf", strata = "population_map_sharks.txt",
  output = c("genepop"))
```

Yet some users may encounter difficulties to install `radiator` package.
In this case, you can use the `vcfR` package instead.
```{r}
vcf_sharks <- read.vcfR("4991_515ind.vcf")
genlight_sharks <- vcfR2genlight(vcf_sharks)
```

Or you can even use the `adegenet` package if you already have a genepop file.
```{r}
genepop_sharks <- read.genepop("4991_515ind.gen", ncode=3L)
```

Download population map file containing the following information: ID, LATITUDE, LONGITUDE and SAMPLING LOCATION.
```{r}
geo <- read.table("population-map-sharks.txt", header=TRUE, sep="\t")
```

## 2 - Find the number of clusters in your dataset

Find clusters into your data usin a genlight or genepop object.
```{r}
grp <- find.clusters(genepop_shark, max.n.clust=10)
```

Choose the number of PCs to retain in relation to the number of individuals in the dataset. This number do not exceed N/3. Keep in mind that more PCs you select, more able you will be to detect clustering in the data.

Choose the best number of K (the smallest BIC value) to keep according to the minimum Cross-validation errror rate gave by the Bayesian Information Criterion. 
Here, a K equal 3 seems to be the best K to keep.

Check the group size of each genetic cluster found.
The step will give you an idea of the releavance of each genetic cluster found.
```{r}
grp$size
```

Save and export the results in a dataframe.
```{r}
data_kmeans <- data.frame(pop$STRATA,grp$grp)
names(data_kmeans)
colnames(data_kmeans) = c("SITE","K")
write.table(data_kmeans, "Individuals_clusters_BIC_515ind.txt", quote=F)
```

[]!(BIC_515ind.png)

Match geographic coordinates of each individual to its inferred genetic group in order to see if there is any clustering relatively to the sampling location.
```{r}
kmean_geo <- merge(data_kmeans, geo, by="IND")
```

Save the results
```{r}
write.table(individuals, "Individuals_clusters_BIC_serranus.txt", quote=F)
```

## 3 - Represent the genetic clusters observed on a sampling map

Download the background data required for creating a map.
```{r}
wH <- map_data("worldHires", xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea
```

Make a map using `ggplot` and `sf` package.
```{r}
x_title="Longitude"
y_title="Latitude"
graph2 <- ggplot() +
  geom_polygon(data = wH, aes(x=long, y = lat, group = group), fill = "gray80", color = NA) +
  coord_fixed(xlim = c(-6, 8), ylim = c(35, 45), ratio=1.2)+
  # facet_wrap(~MEM)+
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(aes(x = LON, y = LAT,fill=K), data=kmean_geo,size=2, shape=21)+
  theme_bw()+theme(legend.position = "none",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y=y_title)+  
  labs(x=x_title)+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))+
  scale_fill_viridis(discrete=TRUE)
#  scale_shape_manual(values=c(21,24))
graph2
```

Save the graph.
```{r}
ggsave("Kmean_map_sharks.pdf")
```

## 4 - Run a DPCA without any prior (based on the K clustering found by the BIC analysis)

Run the DAPC with prior using the number of PCs equal to N/3.
```{r}
dapc_noprior <-dapc(data, grp$grp)
```

Visualize quickly the DPCA results.
```{r}
scatter(dapc2, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
```

Explore the contribution of each SNP to each axis sans save this information.
```{r}
names(dapc2)
head(dapc2$var.contr)
load_dpca2 <- as.data.frame(dapc2$var.contr)
write.table(load_dpca2, "LOADING_prior.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

Observe how much percent of genetic variance is explained by each axis and save this information.
```{r}
pdf("Percent_dapc_all_loci.pdf")
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))
dev.off()
percent
```

## 4 - Visualize the DPCA results with the ggplot package

Save the results of the DPCA in a dataframe.
```{r}
tab=as.data.frame(dapc2$ind.coord)
write.table(tab, "DCPA_results_diplodus.txt", quote=F, sep="\t", row.names=TRUE)
```

Add information to the tab results of the DPCA.
```{r}
tab$IND <- row.names(tab)
```

Add geographical information.
```{r}
dpca_geo <- merge(x = tab, y=geo, by=c("IND"))
```

Make a ggplot graph representing the DAPC for the first and second axes for the regions

```{r}
g = ggplot(dpca_geo, aes(x=LD1, y=LD2, fill=K))+ geom_point(size=2, pch=21)+
  scale_fill_viridis(discrete=TRUE)+
  #guides(fill=FALSE)+
  labs(x="DPC1 (22.6%)")+
  labs(y="DPC2 (18.2%)")+
  theme(axis.text.x=element_text(colour="black"))+
  theme(legend.title=element_blank())+
  theme(axis.text.y=element_text(colour="black",size=12))+
  theme(axis.text.x=element_text(colour="black",size=12))+
  theme(panel.border = element_rect(colour="black", fill=NA, size=3),
        axis.title=element_text(size=18,colour="black",family="Helvetica",face="bold"))+
  theme_bw()+
  theme(legend.title=element_blank())
g
```

Save the ggplot graph.
```{r}
ggsave("DPCA_noprior.pdf",width=13,height=9,dpi=600,units="cm",useDingbats=F)
```

## 5 - Run a DPCA with prior (based on your sampling locations

Find optimal alpha value, which represents the trade-off between power of discrimination and over-fitting. When overfitting, you can observe some clustering even when it is random groups.

```{r}
dapc_a_score <- dapc(data,n.da=156,n.pca=156)
temp_score <- optim.a.score(dapc_a_score)
names(temp_score)
```

Run the DPCA with the number of PCs indicated by the alpha score. 
If you retain too many PCs, you will overfit your data so this step is really crucial
```{r}
dapc_prior <-dapc(data, data@pop)
```
