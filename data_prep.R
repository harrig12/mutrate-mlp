library(rtracklayer, quietly = T)
library(maftools, quietly = T)
library(TCGAmutations, quietly = T)
library(viridis, quietly = T)
library(ggplot2, quietly = T)
library(dplyr, quietly = T)
np <- reticulate::import("numpy")

WINDOW <- 1e6
PEAK_THRESH <- 2

# tensorsigs constants
load(url("http://193.62.55.163/file/R/constants.RData"))
seqlevelsStyle(TS) <- "UCSC"
seqlevelsStyle(RT) <- "UCSC"
seqlevelsStyle(EPI) <- "UCSC"
seqlevelsStyle(NUC) <- "UCSC"

# helpers
get_mafRanges <- function(maf){
    mafRanges <- maf@data %>%
        arrange(Chromosome, Start_Position, End_Position) %>%
        makeGRangesFromDataFrame(keep.extra.columns=T,
                                 ignore.strand=T,
                                 seqinfo=NULL,
                                 seqnames.field=c("Chromosome"),
                                 start.field="Start_Position",
                                 end.field=c("End_Position"),
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)
    seqlevelsStyle(mafRanges) <- "UCSC"
    seqlengths(mafRanges) <- seqlengths(TS)
    return(mafRanges)
}

get_counts <- function(data, tiles, new_col_name = 'count', as_granges = F){
    hits <- findOverlaps(tiles, data)
    counts <- hits %>% 
              as.data.frame() %>% 
              group_by(queryHits) %>% 
              summarize(subjectHits = n()) 

    counts <- counts$subjectHits[match(1:length(tiles), counts$queryHits)]
    counts <- ifelse(is.na(counts), 0, counts)
    if (as_granges){
        mcols(tiles)[[new_col_name]] <- counts
        return(tiles)
    }
    else{ return(counts) }
}


# retrieve bigwigs
tiles <- unlist(tileGenome(seqlengths(TS), tilewidth = WINDOW))
p <- "https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/"
bigwigs <- read.table('bigwigs.txt')$V1
#sel <- grep(pattern = 'E075|E076|E0106|E092|E110|E111', x = bigwigs)
#bigwigs[sel]

epi <- paste0(p, bigwigs)

hist_counts <- list()

for (datasource in epi){
    feature_gr <- rtracklayer::import(BigWigFile(datasource), selection = tiles)
    hist_counts[[datasource]] <- get_counts(feature_gr[mcols(feature_gr)$score > PEAK_THRESH], tiles)
    print('retrieved bigwig')
}

# retrieve mafs
targets = list()
studies = c('COAD', 'STAD')

for (s in studies){
    maf <- TCGAmutations::tcga_load(study = s)
    mafRanges <- get_mafRanges(maf)
    targets[[s]] <- get_counts(mafRanges, tiles)
    print('retrieved maf')
}

# save 

# remove names so that not saved as dict
names(hist_counts) <- names(targets) <- NULL
np$savez('data_all_w1e6.npz', hist_counts=hist_counts, targets=targets)