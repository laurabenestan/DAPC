# Discriminant Analysis in Principal Components (DPCA)

This tutorial provides recommendations for applying the Discriminant Analysis of Principal Components (DAPC) to a wide genomic dataset such as RAD-seq or DarT. 

DPCA acts like a Principal Component Analysis (PCA) and will give more importance to intergroup differences rather than the intragroup differenecs. Groups can be specify at the beginning of the analysis in a blind way (hereafter names DPCA with no prior) or based on sampling information (DPCA with prior). This analysis is particularly relevant in a case of marine species, where low genetic differentaition is often found and hamper the possibility to reveal substle but significant genetic differences. The DPCA was developped by Jombart and colleagues through the `adegenet` package. See the [excellent tutorial](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf) done by Thibault Jombart and Cattlin Collins

## 0 - Prepare your R environment

- [R Version 3.6.0](https://cran.r-project.org/)
	* R packages: adegenet, ade4, ggplot2, reshape, stringr, radiator
  
## 1 - Download genomic and spatial datasets

Transform the vcf to a genepop file.
```{r}
snps_genlight <- genomic_converter(data = "../../00-Data/4991_515ind.vcf", strata = "../../00-Data/population_map_sharks.txt",
  output = c("genepop"))
```

Download genomic information.
```{r}
genpop_shark <- read.genepop("../../00-Data/4991_515ind.gen", ncode=3L)
```

Download geographic information
```{r}
geo <- read.table("../../00-Data/distances_to_mpa_no_coteblue.txt", header=TRUE, sep="\t")
geo$IND <- gsub("_","-", geo$IND)
```

## 2 - Find the number of clusters in your dataset

Find clusters into your data.
```{r}
grp <- find.clusters(genpop_shark, max.n.clust=10)
```

Choose the number of PCs to retain in relation to the number of individuals in the dataset. Choose the best number of K (the smallest BIC value) to keep according to the minimum Cross-validation errror rate. 
Here, a K of 2 seems to be the best K to keep.

Check the group size of each genetic cluster found.
```{r}
grp$size
```

Save the result.
```{r}
data_kmeans <- data.frame(pop$STRATA,grp$grp)
names(data_kmeans)
colnames(data_kmeans) = c("SITE","K")
write.table(data_kmeans, "Individuals_clusters_BIC_515ind.txt", quote=F)
```

Match geographic coordinates with K genetic group info obtained with BIC
```{r}
kmean_geo <- merge(data_kmeans, geo, by="IND")
```

Save the results
```{r}
write.table(individuals, "Individuals_clusters_BIC_serranus.txt", quote=F)
```

## 3 - Represent the genetic clusters obserevd in a map

Download the data for creating a map.
```{r}
wH <- map_data("worldHires", xlim=c(-8,37), ylim=c(29.5,47)) # subset polygons surrounding med sea
```

Make a map using `ggplot`package.
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
ggsave("Kmeans_serranus.pdf")
```

## 4 - Run a DPCA with MPA information (prior)

Find optimal alpha value
```{r}
dapc_a_score <- dapc(data,n.da=156,n.pca=156)
temp_score <- optim.a.score(dapc_a_score)
names(temp_score)
```

Run the DAPC with prior using the number of PCs retained following the alpha score
```{r}
dapc2 <-dapc(data, grp$grp)
```

Visualize quickly the DPCA results.
```{r}
scatter(dapc2, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)
```

Save the results of the DPCA without a prior in a file
```{r}
names(dapc2)
head(dapc2$var.contr)
load_dpca2 <- as.data.frame(dapc2$var.contr)
write.table(load_dpca2, "LOADING_MPAprior_diplodus.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

Analyse how much percent of genetic variance is explained by each axe
```{r}
pdf("Percent_dapc_all_loci.pdf")
percent= dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab="Percent of genetic variance explained by eigenvectors", names.arg=round(percent,2))
dev.off()
percent
```

## 4 - Visualize the DPCA results with the ggplot package

Write the results of the DPCA.
```{r}
tab=as.data.frame(dapc2$ind.coord)
write.table(tab, "DCPA_results_diplodus.txt", quote=F, sep="\t", row.names=TRUE)
```

Add information to the tab results of the DPCA
```{r}
tab$IND <- row.names(tab)
```

Add geographical information
```{r}
dpca_geo <- merge(x = tab, y=geo, by=c("IND"))
```

Make a ggplot graph representing the DAPC for the first and second axes for the regions

```{r}
g = ggplot(dpca_geo, aes(x=LD1, y=LD2, fill=MPA.NAME_EN))+ geom_point(size=2, pch=21)+
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
ggsave("DPCA_MPA_prior.pdf",width=13,height=9,dpi=600,units="cm",useDingbats=F)
```


