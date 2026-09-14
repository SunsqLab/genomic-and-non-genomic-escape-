# Retain genes altered in at least min.num.frequent samples.
FilterGenomicProfiler <- function(genomic.profiler, min.num.frequent=1)
{
	tmp = rowSums(genomic.profiler)
	
	genomic.profiler = genomic.profiler[tmp >= min.num.frequent, ]

return(genomic.profiler)
}
