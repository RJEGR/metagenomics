# Calculate hamming distance over the query sequence
evalDist <- Vectorize(function(query, ref, ...) {
  mmi <- dada2:::nweval(query, ref, ...) # match, mismatch, indel
  rval <- mmi[[2]] + mmi[[3]] # Usual hamming dist
  if(sum(mmi) < nchar(query)) { # query being ends-gapped out
    rval <- rval + nchar(query) - sum(mmi)
    # Add the query ends-gap overlap
  }
  return(rval)
})


# blast parsing
# and load functions

isHit100 <- function(clust, fn) {
  bb <- read.table(fn, comment.char="#", col.names=c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score"))
  bbHit100 <- bb[bb$identity == 100 & bb$coverage == nchar(clust[match(bb$seqid,clust$id),"sequence"]),]
  return(clust$id %in% bbHit100$seqid)
}

isOneOff <- function(clust, fn) {
  hit <- isHit100(clust, fn)
  bball <- read.table(fn, comment.char="#", col.names=c("seqid", "subject", "identity", "coverage", "mismatches", "gaps", "seq_start", "seq_end", "sub_start", "sub_end", "e", "score"))
  bb <- bball[bball$coverage == nchar(clust[match(bball$seqid,clust$id),"sequence"]),] # Only full length hits
  tab <- tapply(bb$identity, bb$seqid, max)
  tab <- tab[match(clust$id, names(tab))]
  seqlens <- nchar(clust$sequence)
  oneOff <- tab<100 & (abs(tab - 100.0*(seqlens-1)/seqlens) < 0.01)
  oneOff[is.na(oneOff)] <- FALSE # happens if no hits were full coverage
  names(oneOff) <- clust$id # Also drop the name to NA so fix here
  # Also get coverage=length-1 matches
  bb <- bball[bball$coverage == nchar(clust[match(bball$seqid,clust$id),"sequence"])-1,] # Full length-1 hits
  bb <- bb[bb$identity==100,]
  oneOff <- oneOff | clust$id %in% bb$seqid
  # Make sure not a hit
  oneOff[hit] <- FALSE
  return(oneOff)
}

#Calculate the hamming distance to the nearest-larger-neighbor (i.e. more abundant output sequence):

get_ham_nln <- function(unqs.in, band=16) {
  require(dada2)
  unqs.in <- getUniques(unqs.in)
  unqs <- sort(unqs.in, decreasing=TRUE)
  sqs <- names(unqs)
  rval <- sapply(seq(2,length(sqs)), function(i) min(nwhamming(sqs[[i]], sqs[1:i-1], band=band) ))
  rval <- c(NA, rval)
  rval[match(names(unqs.in), names(unqs))] # Back to input ordering
}

#


getAccuracy <- function(df) {
  acc <- rep(is.character(NA), nrow(df))
  acc[df$hit] <- "Exact"
  acc[df$oo] <- "One Off"
  acc[!df$oo & !df$hit] <- "Other"
  acc[df$refdist==0] <- "Reference"
  acc[df$refdist==1] <- "One Off"
  acc <- factor(acc, levels=c("Reference", "Exact", "One Off", "Other"))
  return(acc)
}


