## ---------------------------
##
## Script name: Preprocess gene expression data (batch correct & merging)
##
## Author: Mr. Ashleigh C. Myall
##
## Date Created: 01-11-2020
##
## Copyright (c) Ashleigh C. Myall 2020
## Email: a.myall19@imperial.ac.uk
##
## ---------------------------
##
## Notes: Load in data (list of studies) batch correct, merge, and test batch 
## correction with PCA and gene overlaps. For detailed methods see:
##
## ---------------------------

################################################################################

### Load Libs

################################################################################


library(data.table)
library(sva)
library(hues)
library(dplyr)
library(memoise)
library(lmerTest)

################################################################################

###
###                               Functions
###

################################################################################



#-------------------------------------------------------------------------------

## Combat batch correction function

combat_func <- function(data, batch_class, response_class) {
  if(!missing(response_class)){
    corrected = ComBat(dat = t(data), batch=batch_class, 
                       mod=model.matrix(~1+resp, 
                                        data = data.frame(resp=response_class)), 
                       par.prior = T)
  }else{
    corrected = ComBat(dat = t(data), batch=batch_class, par.prior = T)
  }
  as.data.frame(t(corrected))
}

#-------------------------------------------------------------------------------

## Batch correction GSE

batch_gse = function(x){
  print("Running batch correction over GSE")
  data = as.matrix(x$intensity[,2:(ncol(x$intensity)-1)])
  row.names(data) = x$intensity$ref
  gse = x$ref$GSE
  sample.type = x$intensity$sample.type
  if(length(unique(gse))==1){ # no batch correction for single studies
    bc = data
  }else{
    bc = combat_func(data,gse, sample.type)
  }
  return(list("intensity" = bc, "ref" = x$ref,"tag" = x$tag))
}

#-------------------------------------------------------------------------------

## Batch correction GPL

batch_gpl = function(x){
  print("Running batch correction over GPL")
  data = data.matrix(x$intensity[,-1])
  row.names(data) = x$intensity$ref
  gpl = x$ref$GPL
  sample.type = x$ref$sample.type
  if(length(unique(gpl))==1){ # no batch correction for single platforms
    bc = data
  }else{
    bc = combat_func(data,gpl, sample.type)
  }
  return(list("intensity" = bc, "ref" = x$ref,"tag" = x$tag))
}

#-------------------------------------------------------------------------------

## Merge together studies function

mergeInten = function(dt1,dt2,dt3){
  dt1 = as.data.table(dt1,keep.rownames="ref")
  dt2 = as.data.table(dt2,keep.rownames="ref")
  dt3 = as.data.table(dt3,keep.rownames="ref")

  common_cols = intersect(intersect(colnames(dt1), colnames(dt2)),colnames(dt3))
  
  dt1.int = dt1[, common_cols, with=FALSE]
  dt2.int = dt2[, common_cols, with=FALSE]
  dt3.int = dt3[, common_cols, with=FALSE]
  combined.dt = rbindlist(list(dt1.int,dt2.int,dt3.int))

  return(combined.dt)
}

#-------------------------------------------------------------------------------

## Merge wrapper function

mergeStudies = function(data.list){
  return(list("intensity" = mergeInten(data.list$platform1$intensity,
                                data.list$platform2$intensity,
                                data.list$platform3$intensity),
       "ref" = rbindlist(list(data.list$platform1$ref,
                              data.list$platform2$ref,
                              data.list$platform3$ref)),
       "tag" = "merged_data"))
}


#-------------------------------------------------------------------------------

## Plot PCA function

samplesGroupedPCAPlot = function(data, groupingFactors, xlab = 'PC1', ylab = 'PC2', title) {
  uniqueGroups = unique(groupingFactors)
  colors = iwanthue(length(uniqueGroups))
  plot(data$x,col =colors,pch=20,xlab = xlab, ylab = ylab, main = title)
  legend('topleft', legend = uniqueGroups, col = colors, 
         pch = 20, ncol = 2)
}

#-------------------------------------------------------------------------------

## PCA Wrapper function

# Take two lists containing before and after batch correction intensities and compute
# and visualise first two PCA components.

pca_batch = function(x1,x2){
  x1$intensity = as.data.table(x1$intensity)[, !c("ref","sample.type"), with=FALSE] 
  x2$intensity = as.data.table(x2$intensity)[, !c("ref","sample.type"), with=FALSE] 
  
  print("Computing first PCA")
  pca.x1 = prcomp(data.matrix(x1$intensity))
  print("Computing second PCA")
  pca.x2 = prcomp(data.matrix(x2$intensity))
  
  par(mfrow=c(1,2))
  samplesGroupedPCAPlot(pca.x1, as.factor(bc_2$ref$GPL), 
                        xlab = paste("PC1 (", format(pca.x1$sdev[1], digits = 3), '%)', sep=''), 
                        ylab = paste("PC2 (", format(pca.x1$sdev[2], digits = 3), ')%', sep = ''), 
                        title =  'Before Batch: Coloured by Platform')
  samplesGroupedPCAPlot(pca.x2, as.factor(bc_2$ref$GPL), 
                        xlab = paste("PC1 (", format(pca.x2$sdev[1], digits = 3), '%)', sep=''), 
                        ylab = paste("PC2 (", format(pca.x2$sdev[2], digits = 3), ')%', sep = ''), 
                        title =  'After Batch: Coloured by Platform')
  print("Finished PCA")
  return(list("pca_1"=pca.x1,"pca_2"=pca.x2))
}

#-------------------------------------------------------------------------------

## Function to run anovas

runAnovas <- function(lmf,genes,rdata,retAOV=F,retCol="Pr(>F)"){
  lmf <- update(lmf,x~.)
  res <- apply(genes,1,function(x){
    rds <- data.frame(x=x,rdata)
    aovs <- aov(lmf,data=rds)
    if(retAOV){
      return(aovs)
    }
    aovs.s <- summary(aovs)
    if(length(aovs.s) == 1){
      tmp <- aovs.s[[1]][[retCol]]
      names(tmp) <- rownames(aovs.s[[1]])
    }else{
      tmp <- aovs.s[["Error: Within"]][[1]][[retCol]]
      names(tmp) <- rownames(aovs.s[["Error: Within"]][[1]])
    }
    
    return(tmp)
  })
  return(res)
}

runLmer <- function(lmf,genes,rdata,retLmer=F,retCol="Pr(>F)"){
  require(lmerTest)
  lmf <- update(lmf,x~.)
  res <- apply(genes,1,function(x){
    rds <- data.frame(x=x,rdata)
    aovs <- anova(lmer(lmf,data=rds))
    if(retLmer){
      return(aovs)
    }
    tmp <- aovs[,retCol]
    names(tmp) <- rownames(aovs)
    return(tmp)
  })
  return(res)
}

runAnovas.summary <- function(x,error="Error: Within"){
  return(sapply(x,function(z){
    zs <- summary(z)
    if(length(z) == 1){
      tmp <- zs[[1]][[5]]
      names(tmp) <- rownames(zs[[1]])
    }else{
      tmp <- zs[[error]][[1]][[5]]
      names(tmp) <- rownames(zs[[error]][[1]])
    }
    return(tmp)
  }))
}

p.adjust.anova <- function(x,method="BH",byrow=F){
  if(byrow){
    x.adj <- t(apply(x,1,p.adjust,method))
  }else{
    if(any(rownames(x) == "Residuals")){
      x <- x[-match("Residuals",rownames(x)),]
    }
    dx <- c()
    for(i in 1:nrow(x)){
      dx <- c(dx,x[i,])
    }
    x.adj <- p.adjust(dx,method)
    x.adj <- matrix(x.adj,ncol=ncol(x),nrow=nrow(x),byrow=T)
    rownames(x.adj) <- rownames(x)
    colnames(x.adj) <- colnames(x)
  }
  return(x.adj)
}

getSignificantGenes.anova <- function(x,pval,retNames=F,retLength=T){
  if(retLength){
    gs <- apply(x,1,function(y) length(which(y <= pval)))
    return(gs)
  }
  if(retNames){
    gs <- apply(x,1,function(y) names(which(y <= pval)))
  }else{
    gs <- apply(x,1,function(y) which(y <= pval))
  }
  return(gs)
}	


#-------------------------------------------------------------------------------

## Wrapper function to run fisher & hypergeometric test for gene overlaps


test_batch_genes = function(intens_1,intens_2,ref_1,ref_2){
  # If only one study return(Null)
  if(length(unique(ref_2$GSE))==1){
    print("No batch correction ran - so no test needed")
    return(NULL)
  }
  
  # Anova run function
  anova_run <- function(project){
    if(length(unique(classFactor[rownames(project)])) > 1){
      sub_anova <- runAnovas(x~cls,t(project),data.frame(cls=classFactor[rownames(project)]))
      sub_pvals <- p.adjust.anova(sub_anova,"BH")
      genes <- getSignificantGenes.anova(sub_pvals,0.1,retLength=F,retNames=T)
      if(length(genes) == 0){
        genes <- list(cls=c(),Residuals=c())
      }
    }else{
      print(paste(unique(studyFactor[rownames(project)]),"only has 1 class",unique(classFactor[rownames(project)])))
      genes <- list(cls=c(),Residuals=c())
    }
    return(genes)
  }
  
  # Memorise function to save re-computing time
  by_m <- memoise(by)
  
  # Set and name factors
  studyFactor <- ref_1$GSE
  classFactor <- ref_1$sample.type
  names(classFactor) <- ref_1$ref
  names(studyFactor) <- ref_1$ref
  
  # Convert datatables to data.frames
  if(is.data.table(intens_1)){
    intens_1 = as.data.frame(intens_1)
  }
  
  if(is.data.table(intens_2)){
    intens_2 = as.data.frame(intens_2)
  }
  
  # If row names missing and in dataframe add as row names
  if("ref" %in% colnames(intens_1)){
    row.names(intens_1) = intens_1$ref
  }
  
  if("ref" %in% colnames(intens_2)){
    row.names(intens_2) = intens_2$ref
  }
  
  # Remove reference or sample type columns if in data
  intens_1 = intens_1[ , !(names(intens_1) %in% c("ref","sample.type"))]
  intens_2 = intens_2[ , !(names(intens_2) %in% c("ref","sample.type"))]
  
  
  # Siggens for pre-batch
  print("Analysing pre batch corrected data (For large datatsets this can take some time)")
  
  siggenes_prebatch <- by_m(data.matrix(intens_1),studyFactor,anova_run)
  
  # Siggens for post-batch
  print("Analysing post-batch corrected data (For large datatsets this can take some time)")
  siggenes_postbatch <- by_m(data.matrix(intens_2),studyFactor,anova_run)
  
  
  # Compute gene overlaps
  siggens_overlaps <- lapply(1:length(siggenes_prebatch), function(x){
    int <- length(intersect(siggenes_prebatch[[x]]$cls,siggenes_postbatch[[x]]$cls))
    preb <- length(setdiff(siggenes_prebatch[[x]]$cls,siggenes_postbatch[[x]]$cls))
    posb <- length(setdiff(siggenes_postbatch[[x]]$cls,siggenes_prebatch[[x]]$cls))
    all <- length(union(setdiff(colnames(intens_2),siggenes_prebatch[[x]]$cls),setdiff(colnames(intens_2),siggenes_postbatch[[x]]$cls)))
    
    return(c("intersect" = int, "preb" = preb, "posb" = posb, "all" = all))
  })
  names(siggens_overlaps) <- names(siggenes_prebatch)
  
  
  # Do fisher test to check for overlap
  doFisher <- function(list1,list2,total,verbose=F,list1.name="List1",list2.name="List2",all.int=T){
    int <- length(intersect(list1,list2))
    preb <- length(setdiff(list1,list2))
    posb <- length(setdiff(list2,list1))
    if(all.int){
      all <- length(intersect(setdiff(total,list1),setdiff(total,list2)))
    }else{
      all <- length(union(setdiff(total,list1),setdiff(total,list2)))
    }
    mm <- matrix(c(int,preb,posb,all), nrow=2, ncol=2, byrow=T)
    rownames(mm) <- paste(list1.name,c("Sig","NonSig"))
    colnames(mm) <- paste(list2.name,c("Sig","NonSig"))
    if(verbose)
      print(mm)
    return(fisher.test(mm))
  }
  
  # Do hypergeometric test to check for overlap
  doHyperGeometric <- function(list1,list2,total){
    llist1 <- length(list1)
    llist2 <- length(list2)
    ltotal <- length(total)
    int <- length(intersect(list1,list2))
    pm <- 1-phyper(int,llist1,ltotal-llist1,llist2)
    return(pm)
  }
  
  res <- sapply(1:length(siggenes_prebatch), function(x){
    fisher <- doFisher(siggenes_prebatch[[x]]$cls,siggenes_postbatch[[x]]$cls,colnames(intens_2),T,list1.name="Pre-Batch",list2.name="Post-Batch")
    hyper <- doHyperGeometric(siggenes_prebatch[[x]]$cls,siggenes_postbatch[[x]]$cls,colnames(intens_2))
    return(c(fisher=fisher$p.value,hyper=hyper))
  })
  
  
  
  # Plot output
  resDf <- as.data.frame(res)
  colnames(resDf) <- names(siggenes_postbatch)
  resDf <- t(resDf)
  
  par(mfrow=c(2,1))
  barplot(-log10(resDf[,1]+(1e-16)),main = "Fisher Test p values",las=2)
  barplot(-log10(resDf[,2]+(1e-50)),main = "Hyper Geometric Test",las=2)
  
  # Return results
  return(list("res" = resDf,"siggenes_prebatch"=siggenes_prebatch,
              "siggenes_postbatch"=siggenes_postbatch,
              "siggens_overlaps"=siggens_overlaps))
  
}


################################################################################

##########
##########                       Work Flow 
##########

################################################################################

# ------------------------------------------------------------------------------


## 0. Load in example data (non batch corrected). 

# A list containing 3 platform datasets from GPL570,GPL571,GPL9188
data.list = readRDS("./example data/example_raw.RDS")

## 1. Run batch correction over list of studies 
data.list.bc = lapply(data.list, batch_gse)


## 2. Merge batch corrected studies over common columns
bc.merged = mergeStudies(data.list.bc)


## 3. Run batch correction over combined platform
bc_2 = batch_gpl(bc.merged)


## 4. Batch correction tests

# 4.1 - PCA before/after batch correction
pca.anal = pca_batch(bc.merged,bc_2)
#pca.anal2 = pca_batch(data.list$platform1,data.list.bc$platform1)


# 4.2 Test gene overlaps before/after batch correction

# Example test of batch correction 1
gene_test_b1 = test_batch_genes(data.list$platform1$intensity,
                                data.list.bc$platform1$intensity,
                                data.list$platform1$ref,
                                data.list.bc$platform1$ref)

# Example test of batch correction 2
gene_test_b2 = test_batch_genes(bc.merged$intensity,
                                bc_2$intensity,
                                bc.merged$ref,
                                bc_2$ref)

# 5. Save example batch corrected data
saveRDS(bc_2,"./Output/bc_example.RDS")
