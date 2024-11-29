##Citation
#exportGroupBW RCode adopted and modified from Bastien Herve (Baboon61) 
#Apr 4, 2023. https://github.com/stuart-lab/signac/issues/28
###

library("ggplot2")
library("kableExtra")
library("knitr")
library("htmltools")
library("tidyverse")
library('future.apply')
library("dplyr")
library("Seurat")
library('magrittr')
library("Signac")
library('ArchR')

option_list <- list(
  make_option(c("-i", "--indir"), type="character", default = getwd(), help="Directory of seurat Object to read in"),
  make_option(c("-o", "--outdir"), type="character", default = getwd(), help="Directory to Output Group Files"),
  make_option(c("-t", "--tmpdir"), type="character", default = getwd(), help="Directory to store sub_fragments"),
  make_option(c("-n", "--NormMethod"), type="character", default = "TSS.enrichment", help="'TSS.enrichment','ncells','none' or any quantitative values from @meta.data"),
  make_option(c("-a", "--assay"), type="character", default ="ATAC", help="Assay used for fragment files"),
  make_option(c("-g","--genome"),type ="character", default = "hg38", help = "Reference Genome to align reads"),
  make_option(c("-b","--groupBy"),type ="character", default = "hg38", help = "Metadata used to split fragments"),
  make_option(c("-m","--minCells"),type ="integer", default = "5", help = "MinCells Group must have"),
  make_option(c("-c","--chromosome"),type ="character", default = "primary", help = "primary or all"),
  make_option(c("-T","--tileSize"),type ="integer", default = 25, help = "Size of the BW tiles plotted")
  
)

args <- parse_args(OptionParser(option_list=option_list))

if (args$indir == getwd()){
  var = readline(prompt = paste("Is this OBJECT location correct? (Y/N)",args$indir));
  if (var == N){
    stop("Enter a different directory", call.=FALSE)
  }
}if(args$outdir == getwd()){
  var = readline(prompt = paste("Is this OUTPUT location correct? (Y/N)",args$outdir));
  if (var == N){
    stop("Enter a different directory", call.=FALSE)
  }
}if(args$tmpdir == getwd()){
  var = readline(prompt = paste("Is this TMP location correct? (Y/N)",args$tmpdir));
  if (var == N){
    stop("Enter a different directory", call.=FALSE)
  }
}

######Funnction to generate Fragment files
exportGroupBW  <- function(
    object = NULL,
    assay = args$assay,
    genome = args$genome,
    groupBy = args$groupBy,
    downGroupBy = "all",
    normMethod = args$NormMethod, #'ncells', TSS.enrichment ,'none' or any quantitative values from @meta.data
    tileSize = args$tileSize,
    cutoff = 4,
    minCells = args$minCells, #Minimum number of Cells in grouping needed for it to form a fragment file``
    chromosome = args$chromosome, #primary or all
    threads = NULL,
    outdir = args$outdir,
    tmpdir = args$tmpdir,
    group = NULL
){
  
  DefaultAssay(object) <- assay
  genome(object) <- genome
  annot <- object@meta.data
  
  #Column to group the cells by
  Groups <- annot[, groupBy, drop=FALSE]
  cellGroups <- lapply(split(rownames(Groups), factor(apply(Groups, 1, paste, collapse="_"))), unique)
  
  if(!is.null(minCells)){
         cellGroups <- Filter(function(x) length(x) >= minCells, cellGroups)
  }
  cellGroups_rev_list <- unlist(cellGroups)
  names(cellGroups_rev_list) <- rep(names(cellGroups), sapply(cellGroups, length))
  GroupsNames <- unique(names(cellGroups_rev_list))
  
  #GroupsNames <- unique(annot[,groupBy])
  #GroupsNames <- unique(apply(annot[, groupBy], 1, paste, collapse = "_"))
  
  if(unlist(downGroupBy)[1] != 'all'){
    downGroupsNames <- apply(do.call(expand.grid, downGroupBy),1, paste, collapse = "_")
    GroupsNames <- GroupsNames[GroupsNames %in% downGroupsNames]
  }
  
  #Column to normalized by
  if (normMethod == 'ncells'){
    normBy <- normMethod
  } else{
    normBy <- annot[, normMethod, drop=FALSE]
  }
  
  #Get chromosome information
  if(chromosome=="primary"){
    prim_chr <- names(seqlengths(object)[!grepl("_alt|_fix|_random|chrUn", names(seqlengths(object)))])
    seqlevels(object) <- prim_chr
  }
  
  availableChr<- object@assays[["ATAC"]]@annotation@seqinfo@seqnames
  
  ChromLengths<- object@assays[["ATAC"]]@annotation@seqinfo
  chromLengths <- seqlengths(ChromLengths)
  
  chromSizes <- GRanges(seqnames = availableChr, IRanges(start = rep(1, length(availableChr)), end = as.numeric(chromLengths)))
  
  #Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(slidingWindows(chromSizes, width = tileSize, step = tileSize))
  
  if (threads == 1){
    avai_threads <- 1
  }else{
    avai_threads <- threads - 1
  }
  
  #Set number of thread in future
  plan("multicore", workers = avai_threads)
  
  #Run the creation of bigwig for each cellgroups
  covFiles <- future_lapply(GroupsNames, .createBWGroup, availableChr, tiles, normBy, tileSize, normMethod, cutoff, outdir, tmpdir)
  
  plan("multicore", workers = 1)
  
  covFiles
  
}

.createBWGroup <- function(groupNamei, availableChr, tiles, normBy, tileSize, normMethod, cutoff, outdir, tmpdir){
  
  #Read the fragments file associated to the group
  fragi <- rtracklayer::import(paste0(tmpdir, "/", groupNamei, ".tsv.gz"),format = "bed")
  
  cellGroupi <- unique(fragi$name)
  
  covFile <- file.path(outdir, paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw"))
  
  covList <- lapply(seq_along(availableChr), function(k){
    
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    
    tilesk <- tiles[BiocGenerics::which(seqnames(tiles) %bcin% availableChr[k])]
    
    if(length(fragik) == 0){
      tilesk$reads <- 0
      #If fragments
    }else{
      
      #N Tiles
      
      nTiles <- chromLengths[availableChr[k]] / tileSize #original
      
      #Add one tile if there is extra bases
      if (length(nTiles)%%1 != 0) { #had to add length
        nTiles <- trunc(nTiles) + 1
      }
      
      #Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      
      #For each tiles of this chromosome, create start tile and end tile row, set the associated counts matching with the fragments
      
      mat <- Matrix::sparseMatrix( 
        i = c(trunc(start(fragik) / tileSize), trunc(end(fragik) / tileSize)) + 1,
        j = as.vector(c(matchID, matchID)),
        x = rep(1,  2*length(fragik)),
        dims = c(tilesk@elementMetadata@nrows, length(cellGroupi))) #Had to change from nTiles for tilesk as it seems sparse matrix becomes its input
      
      #Max count for a cells in a tile is set to ceiling (4)
      if(!is.null(cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      
      #Sums the cells
      mat <- Matrix::rowSums(mat)
      
      tilesk$reads <- mat
      
      #Normalization of counts by the sum of readsintss for each cells in group
      if(normMethod == "ncells"){
        tilesk$reads <- tilesk$reads*10^4 / length(cellGroupi)
      }else if(tolower(normMethod) %in% c("none")){
      }else{
        tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
      }
    }
    
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    
    tilesk
    
  })
  
  names(covList) <- availableChr
  
  covList <- as(covList, "RleList")
  
  rtracklayer::export.bw(object = covList, con = covFile)
  
  covFile
}

atac_seurat<- readRDS(args$indir)

#Run function definitions before running this function
exportGroupBW(atac_seurat, group = group)

