# Discriminant Analysis in Principal Components (DPCA)

This tutorial provides recommendations for applying the Discriminant Analysis of Principal Components (DAPC) to a wide genomic dataset such as RAD-seq or DarT. 

DPCA acts like a Principal Component Analysis (PCA) and will give more importance to intergroup differences rather than the intragroup differenecs. Groups can be specify at the beginning of the analysis in a blind way (hereafter names DPCA with no prior) or based on sampling information (DPCA with prior). This analysis is particularly relevant in a case of marine species, where low genetic differentaition is often found and hamper the possibility to reveal substle but significant genetic differences. The DPCA was developped by Jombart and colleagues through the `adegenet` package. See the [excellent tutorial](http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf) done by Thibault Jombart and Cattlin Collins

## 0 - Prepare your R environment

- [R Version 3.6.0](https://cran.r-project.org/)
	* R packages: adegenet, ade4, ggplot2, reshape, stringr, radiator
  
## 1 - Download dataset
pop <- read.table("../../00-Data/population_map_sharks.txt", header=TRUE, sep="\t")
pop_strata_vector <- levels(pop$STRATA)

Import vcf in genind
```{r}
snps_genlight <- genomic_converter(data = "../../00-Data/4991_515ind.vcf", strata = "../../00-Data/population_map_sharks.txt",
  output = c("genlight"))
```

```{r}
genpop_shark <- read.genepop("../../00-Data/4991_515ind.gen", ncode=3L)
```

Find clusters into your data.
```{r}
grp <- find.clusters(genpop_shark, max.n.clust=10)
```

Choose the number of PC to retain in relation to the number of individuals in the dataset.
Choose the number of K (the smallest BIC value).

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



