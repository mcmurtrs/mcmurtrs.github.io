# Creating a Maximum Likelihood Phylogenetic Tree with Randomized Axelerated Maximum Likelihood (MAxML)

#### Authors: Shawn McMurtrey, Carolina Piña Páez, Nick Carleson, Patrick Bennett, Ricardo I. Alcalá Briseño, Javier Tabima, Jared LeBoldus

## Step 1: Convert VCF file into fasta file (Nexus format) with vcf2phylip python script:
- Before we can run RAxML we need to know which evolutionary model to use.
- This can be figured out with JMODELTEST however in order to run JMODELTEST we need a Nexus or fasta file.
- We only have a vcf file which doesn't contain sequence information. In order to make a P. tree we will need sequence information.

```
#On the command line download the python script by entering the following:
wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py

#Change permissions of the python file so that it is executable
chmod +x vcf2phylip.py

#Run this command in the directory that you have the vcf file and the python file at. -f creates a fasta file. https://github.com/edgardomortiz/vcf2phylip
python vcf2phylip.py -i Cs_ALL_filtered_NEW_11_30.vcf.gz -f
```

## Step 2: Run evolutionary model test with RAxML-ng to determine best fitting model to use:

```
/nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/bin/modeltest-ng-static -i 30_final_Dec21.min4.fasta -d nt
```

## Step 3: After model-test-ng has finished running, Build ML Starting tree with RAxML-ng
```
SGE_Batch -q bpp@symbiosis -c '/nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/bin/raxml-ng --msa 30_final_Dec21.min4.fasta --model TVM --prefix T1 --threads 20' -r T1 -P 20
```

## Step 4: Make 1000 bootstrap trees
```
SGE_Batch -q bpp@symbiosis -c '/nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/bin/raxml-ng --bootstrap --msa T1.raxml.rba --model TVM --prefix T2 --threads 20 --bs-tree 1000' -r T2 -P 20
```

## Step 5: Check the convergence of the boostrapping test 

```
/nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/bin/raxml-ng --bsconverge --bs-trees T2.raxml.bootstraps --prefix Test --threads 2 --bs-cutoff 0.03
```

## Step 6: Finally map Support values from Bootstrap test to best scoring ML tree
- After this step is finished it will print out where our final tree is at!
- i.e. Best ML tree with Felsenstein bootstrap (FBP) support values saved to: /nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/cs_align/Phylo_tree/T3.raxml.support
- Hooray!
```
/nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/bin/raxml-ng --support --tree T1.raxml.bestTree --bs-trees T2.raxml.bootstraps --prefix T3 --threads 2
```

## Step 7: Move files from the cluster to your computer and change file extension
- /nfs1/BPP/LeBoldus_Lab/user_folders/mcmurtrs/cs_align/Phylo_tree/T3.raxml.support
- Move the .raxml.support file to your desktop
- This file contains the tree and bootstrap values.
- On your computer, change the extension of the bipartions file to .tre


## Step 8: Read the .tre file into R

```{r}

library(ggtree)
library(ape)

myTree <- read.tree("C:/Users/shawn/Desktop/RAxML_bipartitions.tre")

```

## Step 9: Open the tree file and inspect the contents

- Notice that this is a new type of file.
- It is not a data file as we would normally expect to see in R, this is a list.
- This particular list file contains 5 elements: edge, edge.length, Nnode, node.label, and tip.label
- For this example, there are 173 edges, 173 edge lengths, 86 nodes, 86 node lables and 88 node tip labels. 

![image](https://user-images.githubusercontent.com/49656044/144702891-56af2f7d-30b4-4dbe-ab1e-b5da1a2f931d.png)


## Step 10: Reroot the tree

```{r, message=FALSE}
library(ape)
library(phangorn)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

myTree <- midpoint(myTree)

```

## Step 11: Add labels with fortify (not sure how to do this yet...)

```{r}

myTree.fort <- fortify(myTree)
myTree.fort$bootstrap <- NA
myTree.fort$bootstrap[!myTree.fort$isTip]<-myTree$node.label

```


## Step 12: Filter bootstrap labels to remove those less than 70

```{r}
q3 <- ggtree(myTree.fort)
d3 <- q3$data
d3 <- d3[!d3$isTip,]
d3$label <- as.numeric(d3$label)
d3 <- d3[d3$label > 70,]
```

## Step 13: I'm not really sure why we do this step??


```{r}
RAXtre<-q3
```



## Step 14: Visualize the tree

```{r}
# The tree
RAXtre + 

# Add ML values > 70
geom_nodelab(aes(y=branch,  label = label, subset =   !is.na(as.numeric(label)) & as.numeric(label) >70), vjust=-.5, size=.1) +
  
# Tree tip labels
geom_tiplab(size = 2) +
  
# Bootstrap values  
geom_nodepoint(data=d3,size=2,shape=21,aes(fill=cut(as.numeric(bootstrap),c(70,80,90,100),left=T, include.lowest=T))) +
  
# Override colors and turn to grey scale 
scale_fill_manual(name = "Bootstrap support",values = c("white","grey60","black"),labels = c("70-79", "80-89", "90-100"))+
  
# Legend 
guides(colour=guide_legend(override.aes=list(size=10,shape=16,alpha=1))) +
  
# Legend position
theme_tree(legend.position=c(0.6,0.75),legend.text=element_text(size=12),legend.title=element_text(size=12),text=element_text(size=12,face="bold"),legend.margin=margin(t = 0, unit='cm'),legend.box.margin=margin(t = 0, unit='cm'))+ 
  
# Tree scale bar
geom_treescale(x=.5,y=-15,width=15,fontsize=3,offset=0.5)

```





#### Helpful githubs/tutorials:

http://treethinkers.org/tutorials/model-selection/

https://knausb.github.io/vcfR_documentation/matrices.html

http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/species-tree-inference/

https://github.com/Patrick-Bennett/leptographium_popgen/blob/master/PhylogeneticTree/RAxML_Tree_LWP_LWW_FINAL.Rmd

https://evomics.org/learning/phylogenetics/jmodeltest/



#### Helpful videos:

#Installing Apache Ant
https://www.youtube.com/watch?v=7z2yXY57jxY

#Info on Phylogenetic Trees
https://www.youtube.com/watch?v=uAf0DCGynTk&t=323s

#JModelTest Info
https://www.youtube.com/watch?v=uOWxBDUcss4&t=366s


 

