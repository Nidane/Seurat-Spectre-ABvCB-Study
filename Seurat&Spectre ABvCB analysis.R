####
#### 0. Analysis session setup #################################################
####

### Load packages
library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

library(tidyverse)
library(dplyr)
library(ggplot2)
library(CytoNorm)
library(flowCore)
library(Seurat)
library(SeuratObject)
library(ggthemes)
library(devtools)
library(stats)

### Set primary directory
dirname(rstudioapi::getActiveDocumentContext()$path)            # Finds the directory where this script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Sets the working directory to where the script is located
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory
setwd("C:/Users/Data")
setwd("CSV/") # Sample files exported from FlowJo as .csv with scaled values into this folder
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory
setwd("C:/Users/Data")
setwd("METADATA/") # Contains sample information and groups
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory
dir.create("OUTPUT FOLDER")
setwd("OUTPUT FOLDER")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

### Import data
setwd(InputDirectory)
list.files(InputDirectory, ".csv")

data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)

# Merge sample data into a single matrix
data.list[[1]]
cell.dat <- Spectre::do.merge.files(dat = data.list)

##### Choose essential data and rename the markers ####
as.matrix(names(cell.dat))
cell.dat = cell.dat[,c(8:27,31:32)]

#### 1. ArcSinh Transformation (using Spectre) #################################
setwd(OutputDirectory)
dir.create("1 Arcsinh")
setwd("1 Arcsinh")

##### Choose markers to transform with arcsinh #####
to.asinh <- names(cell.dat)[c(1:20)]
to.asinh

cofactor <- 2000
cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)

##### Remove non-transformed columns #####
cell.dat = cell.dat[,c(21:42)]


#### 1-2 Add Metadata ##########################################################
setwd(MetaDirectory)

# Add metadata to data.table
meta.dat <- fread("sample.details.csv")
cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)

### Now that the data is prepared, we can perform;
# 2. general PCA analysis based on the MFI of all samples
# 3. Cluster analysis with Seurat 
# OR 
# 4. Cluster analysis with Spectre

#### 2. PCA based on MFI #######################################################
setwd(OutputDirectory)
dir.create("2 PCA")
setwd("2 PCA")

data = cell.dat

# Choose markers to perform PCA
exp <- data[,c(3:7,9:10,12:21)] 

exp.t <- t(exp)
exp.t.pc <- prcomp(exp, center=TRUE, scale=TRUE)

pc <- exp.t.pc$x
pc.df = as.data.frame(pc)

#### 3. SEURAT ANALYSIS ########################################################
setwd(OutputDirectory)
dir.create("Seurat")
setwd("Seurat")

# Reformatting before creating seurat object
pre.seu <- as.data.frame(cell.dat)

# Re set the column names
as.matrix(names(pre.seu))
for ( col in 1:ncol(pre.seu)){
  colnames(pre.seu)[col] <-  sub("_.*", "", colnames(pre.seu)[col])
}
for (col in 1:ncol(pre.seu)){
  colnames(pre.seu)[col] <-  sub("-", "", colnames(pre.seu)[col])
}
colnames(pre.seu)

# Choose markers to use for clustering
pre.seu = pre.seu[,c(1:2,4:8,10:11,14:18,20)]

##### Downsample per Sample ####
pre.seu = pre.seu %>% group_by(Sample) %>% slice_sample(n=7500)

# Extract metadata
pre.seu.meta = pre.seu[,1:2]

# Create the count matrix for the flow data, with individual cells as columns and markers as rows 
pre.seu.obj = pre.seu[,3:ncol(pre.seu)]
rownames(pre.seu.obj) = 1:nrow(pre.seu.obj)
pre.seu.obj = t(pre.seu.obj)

# Create the seurat object (add the meta data and label project)
seuobj = CreateSeuratObject(counts = pre.seu.obj, project = "CD8", assay = "Flow",
                            meta.data = pre.seu.meta)

##### PCA in Seurat #####
setwd(OutputDirectory)
dir.create("Seurat")
setwd("Seurat")
dir.create("PCA")
setwd("PCA")

dat <- seuobj

dat <- Seurat::FindVariableFeatures(dat, nfeatures = nrow(dat))

# pseudo scaling of the data to proceed the Seurat workflow 
# did not actually change the data
VariableFeaturePlot(dat)
dat <- ScaleData(object = dat, do.scale = FALSE, do.center = FALSE)  

# Perform PCA
dat <- RunPCA(dat, npcs = 50)

# Total dimensions
n.dims = ncol(dat@reductions$pca)


##### Clustering and Dimension Reduction #####
# Choose number of dimensions to include based on PCA
dat = FindNeighbors(dat, dims = 1:8)

# Change number of clusters here with resolution
dat <- FindClusters(dat, resolution = 0.5)

# UMAP and TSNE 
dat <- RunUMAP(dat, dims = 1:8)
dat <- RunTSNE(dat, dims = 1:8)

#### 4. SPECTRE ANALYSIS #######################################################
setwd(OutputDirectory)
dir.create("Spectre")
setwd("Spectre")

# Specify the markers for clustering
# T cell lineage markers are excluded and only functional markers are included 
cluster.cols <- names(cell.dat)[c(3:7,9:10,12:18,20:22)] 

##### Downsample per sample #####
data.frame(table(cell.dat[[sample.col]]))
sub.cl = rep(c(7500),each=10)

rm(cell.sub.cl)
cell.sub.cl <- do.subsample(cell.dat, sub.cl, Sample)

##### FlowSOM Clustering #####
cell.sub.cl <- run.flowsom(cell.sub.cl, cluster.cols, meta.k = 16,
                           ydim=5, xdim=5, # can change SOM map size here
)

##### Dimensional Reduction #####
cell.sub.cl <- run.umap(cell.sub.cl, cluster.cols)