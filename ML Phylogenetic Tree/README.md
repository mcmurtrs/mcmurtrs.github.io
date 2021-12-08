# Creating a Maximum Likelihood Phylogenetic Tree with Randomized Axelerated Maximum Likelihood (MAxML)

#### Authors: Shawn McMurtrey, Carolina Piña Páez, Nick Carleson, Patrick Bennett, Ricardo I. Alcalá Briseño, Javier Tabima, Jared LeBoldus

## Step 1: Import filtered VCF file into R
- Read VCF file into R with vcfR package:

```{r}
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
LRR_VCF <- read.vcfR("C:/Users/shawn/Desktop/Cs_ALL_filtered_NEW_11_30.vcf.gz")
LRR_VCF

```

## Step 2: Convert the very large VCF file into a matrix for computational efficeincy 
- Format missing data and account for mixed phases with pipes into a gt matrix.
- Print out the number of samples and number of SNPs within your file
- Open the file in Notepad ++, save those numbers and paste them in the first line of "gt.phy" that you're going to save in the next two step.
- Do a "Find and Replace" to replace all of the " symbols with a blank space (the " symbol throws an error during the next step when you run MAxML on the cluster.

```{r}

gt <- extract.gt(LRR_VCF, element = "GT")
gt[gt == "0/0"] <- 0
gt[gt == "0|0"] <- 0
gt[gt == "0/1"] <- 1
gt[gt == "0|1"] <- 1
gt[gt == "1|0"] <- 1
gt[gt == "1/1"] <- 2
gt[gt == "1|1"] <- 2
gt[is.na(gt)] <- "?"
gt <- t(gt)
dim(gt)
#102 1341746

```

## Step 3: Write out the matrix into a table and save it to a file called gt.phy

```{r}

gt.pyl <- apply(gt, 1, paste0, collapse="")
write.table(gt.pyl, sep = "\t", col.names = F, file = "gt.phy")

```

## Step 4: Run RAxML on the cluster
- Note that these parameters are specifically for the Center for Quantitative Life Sciences cluster at Oregon State University

```
SGE_Batch -q bpp@symbiosis -c 'mpiexec -n 10 raxmlHPC-MPI -N 50 -n myMLJob -s gt.phy -m MULTICAT -f a -x 12345 -p 12345' -r snp_tree -P 10

```

## Step 6: Move files from the cluster to your computer and change file extension

- Move the .bipartitions file to your desktop
- This file contains the tree and bootstrap values.
- On your computer, change the extension of the bipartions file to .tre


## Step 7: Read the .tre file into R

```{r}

library(ggtree)
library(ape)

myTree <- read.tree("C:/Users/shawn/Desktop/RAxML_bipartitions.tre")

```

## Step 8: Open the tree file and inspect the contents

- Notice that this is a new type of file.
- It is not a data file as we would normally expect to see in R, this is a list.
- This particular list file contains 5 elements: edge, edge.length, Nnode, node.label, and tip.label
- For this example, there are 173 edges, 173 edge lengths, 86 nodes, 86 node lables and 88 node tip labels. 

![image](https://user-images.githubusercontent.com/49656044/144702891-56af2f7d-30b4-4dbe-ab1e-b5da1a2f931d.png)


## Step 9: Reroot the tree

```{r, message=FALSE}
library(ape)
library(phangorn)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

myTree <- midpoint(myTree)

```

## Step 10: Add labels with fortify

```{r}

myTree.fort <- fortify(myTree)
myTree.fort$bootstrap <- NA
myTree.fort$bootstrap[!myTree.fort$isTip]<-myTree$node.label

```


## Step 11: Filter bootstrap labels to remove those less than 70

```{r}
q3 <- ggtree(myTree.fort)
d3 <- q3$data
d3 <- d3[!d3$isTip,]
d3$label <- as.numeric(d3$label)
d3 <- d3[d3$label > 70,]
```

## Step 12: I'm not really sure why we do this step??


```{r}
RAXtre<-q3
```



## Step 13: Visualize the tree

```{r}
# The tree
RAXtre +

# Add ML values > 70  
geom_nodelab(aes(y=branch, label=label, subset =   !is.na(as.numeric(label)) & as.numeric(label) >70), vjust=-.5, size=2) +
  
# Tree tip labels
geom_tiplab(size = 2) +
  
# Tree scale bar
geom_treescale(x=.5,y=-15,width=15,fontsize=3,offset=0.5)
```





#### Helpful githubs/tutorials:
https://knausb.github.io/vcfR_documentation/matrices.html

https://github.com/Patrick-Bennett/leptographium_popgen/blob/master/PhylogeneticTree/RAxML_Tree_LWP_LWW_FINAL.Rmd

#### Helpful videos:
https://www.youtube.com/watch?v=uAf0DCGynTk&t=323s




