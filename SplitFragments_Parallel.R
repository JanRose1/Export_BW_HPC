
###CITATION
#SplitFragments RCode adopted and modified from Bastien Herve (Baboon61) 
#Apr 4, 2023. https://github.com/stuart-lab/signac/issues/28
###


library('rtracklayer')
library("knitr")
library('future.apply')
library("dplyr")
library("Seurat")
library("Signac")
library('ArchR')
library('slurmR')
library('lme4')
library('foreach')
library('doParallel')


atac_seurat<- readRDS("/nfs/turbo/path-rjhryan-turbo/lab-members/Jan/Ryan_Lab_Chromatin/Level_3_ATAC/sc_ATAC_Lv3")


splitGroupFragments <- function(
    object = NULL,
    assay = "ATAC",
    genome = "hg38",
    groupBy = "annotation_20230508", #Metadata used to split fragments
    downGroupBy = "all", #Only change if interested in specific clusters or timepoints. Otherwise use all
    minCells = 5, 
    maxCells = 10000,
    nbFrags = 10000000,
    threads=1, #Keep at 1
    setCores = 1, #Set this value to num of cpus used (-1 if on local machine)
    test=FALSE, #Tests code on a chunk of a fragment file 
    test.nbFrags=NULL,
    outdir = NULL, #Directory to output final fragments
    tmpdir = NULL, #Directory to store sub-fragment files
    frgdir = NULL, #Directory of the Merged fragment file
){
  
  dir.create(file.path(outdir), showWarnings = FALSE)
  dir.create(tmpdir, showWarnings = FALSE)
  
  DefaultAssay(object) <- assay
  genome(object) <- genome
  
  #Get fragments file path
  if (length(Fragments(object)) > 1){
    print("The object contains a list of fragments files")
    #frag_file <- unique(unlist(lapply(Fragments(object), function(x){GetFragmentData(x)})))
    frag_file <- frgdir
    if (length(frag_file) > 1){
      #Could also try to run on each specific fragments files without splitting
      print("The list of fragments files are pointing to different fragment files")
      print("Please, aggregate the fragments files before submitting")
    }else{
      bed_input_path <- frag_file
    }
  }else{
    bed_input_path <- GetFragmentData(Fragments(object))
  }
  
  if(test){
    dir.create(paste0(tmpdir, "/test"), showWarnings = FALSE)
    nbFrags <- test.nbFrags
    print(paste0("Testing on a small chunk of the fragments file, selecting : ",nbFrags, " fragments"))
    system(sprintf("zcat %s | head -%d > %s/test/test_%s_atac_fragments.tsv", bed_input_path, nbFrags, tmpdir, as.character(nbFrags)))
    system(sprintf("gzip %s/test/test_%s_atac_fragments.tsv", tmpdir, as.character(nbFrags)))
    bed_input_path <- paste0(tmpdir, "/test/test_", nbFrags, "_atac_fragments.tsv.gz")
    print(paste0("Testing file located at : ",bed_input_path))
  }
  
  annot <- object@meta.data
  
  #Column to group the cells by
  Groups <- annot[, groupBy, drop=FALSE]
  
  if(unlist(downGroupBy)[1] != 'all'){
    Groups <- Groups[which(Reduce(`&`, Map(`%in%`, Groups[intersect(names(Groups), names(downGroupBy))], downGroupBy))),]
  }
  
  #Split cell names by group
  cellGroups <- lapply(split(rownames(Groups), factor(apply(Groups, 1, paste, collapse="_"))), unique)
  
  #Sample each cell group to maxCells
  if(!is.null(maxCells)){
    cellGroups <- lapply(cellGroups, function(x){if(length(x) <= maxCells){x}else{sample(x, maxCells)}})
  }
  
  #Remove group with less than mixCells
  if(!is.null(minCells)){
    cellGroups <- Filter(function(x) length(x) >= minCells, cellGroups)
  }
  
  #Create a vector of selected cells with the group as name
  cellGroups_rev_list <- unlist(cellGroups)
  names(cellGroups_rev_list) <- rep(names(cellGroups), sapply(cellGroups, length))
  
  ##If there is 1B line in the fragment file
  ##We can split the fragments file by the number of available threads (minus 1 thread not to overload)
  
  if (threads == 1){
    avai_threads <- 1
  }else{
    avai_threads <- threads - 1
  }
  theorical_nbfile_per_thread <- nbFrags/avai_threads
  suitable_nbfile_per_thread <- ceiling(theorical_nbfile_per_thread)
  print(paste0("Each sub-fragments files will contain : ", suitable_nbfile_per_thread," fragments"))
  system(sprintf("gunzip -c %s | split --lines %d --additional-suffix=.tsv - %s/sub_atac_fragments", bed_input_path, suitable_nbfile_per_thread, tmpdir))
  
  #Set number of thread in future
  plan("multicore", workers = avai_threads)
  
  #Function to generate fragments files
  groupfragments <- function(x, cellname, cellgroup, tmpdir) {
    bed_input <- gzfile(paste0(tmpdir,"/",x), "r")
    bed_output <- file(paste0(tmpdir, "/", cellgroup, "_", strsplit(x,"\\.")[[1]][1], ".tsv"), "w")
    
    while(length(line <- readLines(bed_input, n = 1)) > 0) {
      cell_name_bed <- strsplit(line, "\t")[[1]][4]
      if (cell_name_bed %in% cellname) {
        writeLines(line, bed_output)
      }
    }
    close(bed_input)
    close(bed_output)
  }
  
  #Create a cluster
  numCores <- setCores #One less than avail
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  files <- list.files(tmpdir)[list.files(tmpdir) %like% "sub_atac_fragments"] 
  foreach(group = unique(names(cellGroups_rev_list))) %dopar% {
    if (file.exists(file = paste(tmpdir,paste(group,'.tsv.gz')))){
      next
    }
    cells <- cellGroups_rev_list[names(cellGroups_rev_list)==group]
    future.apply::future_lapply(files,groupfragments, cells, group, tmpdir, future.stdout=FALSE)
    system(sprintf("cat %s/%s* | sort -k1,1 -k2,2n -k3,3n - > %s/%s.tsv", tmpdir, group, tmpdir, group))
    system(sprintf("rm %s/%s_sub*", tmpdir, group))
    system(sprintf("gzip %s/%s.tsv", tmpdir, group))
  }
  
  stopCluster(cl)
  
  #Removing sub-fragments files
  system(sprintf("rm %s/sub_atac_fragments*", tmpdir))
  
  plan("multicore", workers = 1)
  
  if(test){
    system(sprintf("rm -r %s/test", tmpdir))
  }
}


splitGroupFragments(atac_seurat, assay = 'ATAC', genome = 'hg38', threads = 1  , outdir = "/nfs/turbo/path-rjhryan-turbo/lab-members/Jan/Ryan_Lab_Chromatin/Re-Fragments_Annot_20230508", tmpdir = "/nfs/turbo/path-rjhryan-turbo/lab-members/Jan/Ryan_Lab_Chromatin/Re-Fragments_Annot_20230508")
