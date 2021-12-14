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

## Step 4: Convert VCF file into fasta file (Nexus format) with vcf2phylip python script:
- Before we can run RAxML we need to know which evolutionary model to use.
- This can be figured out with JMODELTEST however in order to run JMODELTEST we need a Nexus or fasta file.
- We only have a vcf file which doesn't contain sequence information. In order to make a P. tree we need sequence information.

```
#On the command line download the python script by entering the following:
wget https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py

#Change permissions of the python file so that it is executable
chmod +x vcf2phylip.py

#Run this command in the directory that you have the vcf file and the python file at. -f creates a fasta file. https://github.com/edgardomortiz/vcf2phylip
python vcf2phylip.py -i Cs_ALL_filtered_NEW_11_30.vcf.gz -f
```

## Step 5: Run JMODELTEST to determine the best fitting evolutionary model to use:
- In order to install the GUI of JMODELTEST you need to first intall Apache Ant and Java JDK 1.6 (or a newer version).
- This video does an excellent job of showing how to install Apache Ant on Windows 10: https://www.youtube.com/watch?v=7z2yXY57jxY
- One you have correctly installed Apache Ant, Java JDK 1.6, and JMORDELTEST you can then open the GUI of JMODELTEST with by selecting the .jar file and then follow the tutorial found here for running the JMODELTEST: https://evomics.org/learning/phylogenetics/jmodeltest/
- Note that this step is computationally heavy and appears to be taking a long time to run.
- As a trial run, I ran a fasta file containing 20 full genome samples and it took ~18 hours to run.
- The large file that contained 101 samples would not even load into the GUI so an alternative method will need to be identified for using larger files. 

## Step 5.1: After JModelTest has Finished Running...

- The results can be found by watching the video here:
- https://www.youtube.com/watch?v=uOWxBDUcss4&t=366s
- In short, 'Analysis' > 'Do BIC Calculations' > 'Keep default settings unless you have a reason not to' > 'Results' > 'Show results table'


## Analysis of Results

- We can see from the screenshots below that the best fitting evolutionary model from this analysis is the GTR model.
- We know this because it has the highest AIC, BIC, and DT values and also has a value of delta = 0 for all comparisons.

![image](https://user-images.githubusercontent.com/49656044/146057741-db8ef153-9ead-4375-a403-791b0f1b3585.png)

![image](https://user-images.githubusercontent.com/49656044/146057683-6ac36d70-b8cc-4941-a101-05be0942d95e.png)

![image](https://user-images.githubusercontent.com/49656044/146057791-d63909fc-c215-4102-9dd2-913baa49a996.png)


## Full JModelTest Results for 20 Samples from scattered places across North American, Siberia, and Japan 

https://mcmurtrs.github.io/lrr.fastqc.github.io/Final_all_over_dec12.min4.fasta.jmodeltest.html




## Step 6: Run RAxML on the cluster
- Now that we know which evolutionary model to use, we can change this in the one liner command for RAxML and start builing our tree!
- Note that the parameters below are specifically for the Center for Quantitative Life Sciences cluster at Oregon State University
- Instead of just blindly copying and pasting the command below (as I did the first time *rolls eyes*) lets dissect the contents and try to understand what each parameter means.
- 'mpiexec' is a the respective Message Passing Interface (MPI) run-time command, e.g. mpiexec or mpirun depending on your
local installation (please check with your local computer scientist or someone at the CQLS).
- 'raxmlHPC-MPI' the version of RAxML being used; From the manual, RAxML-VI-HPC: Improvement of the OpenMP parallelization, technical
tuning of the GTR-based likelihood models, re-implementation of the MP (Maximum Parsimony) starting tree computations, MPI-based parallel version for multiple bootstrapping.
- '-N' Specifies the number of alternative runs on distinct starting trees, e.g., if '-N 10' is specified, RAxML will compute 10 distinct ML trees starting from 10 distinct randomized maximum parsimony starting trees 
- '-n' specifies the name of the output file
- '-s' is the input the file 
- '-m' is the model to use (Generally, CAT should perform better than GAMMA, particularly with respect to inference times. CAT can typically also produce better scoring trees.)
- '-x' = rapidBootstrapRandomNumberSeed (whatever that means)
- '-p' = parsimonyRandomSeed (again at a loss for what this interprets to)


```
#Example 1, using GTR Model
SGE_Batch -q bpp@symbiosis -c 'mpiexec -n 10 raxmlHPC-MPI -N 50 -n myMLJob -s gt.phy -m GTRCAT -f a -x 12345 -p 12345' -r snp_tree -P 10

#Example 2, using MULTICAT Model 
SGE_Batch -q bpp@symbiosis -c 'mpiexec -n 10 raxmlHPC-MPI -N 50 -n myMLJob -s gt.phy -m MULTICAT -f a -x 12345 -p 12345' -r snp_tree -P 10

```

## Step 7: Move files from the cluster to your computer and change file extension

- Move the .bipartitions file to your desktop
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






