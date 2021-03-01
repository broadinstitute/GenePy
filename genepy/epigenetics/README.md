# epigenomics

Especially targeted to functions related to the analysis of epigenomics data. It has functions to read, merge, denoise, ChIP seq data.

## Available functions:

_in ./chipseq.py_

- bigWigFrom
- ReadRoseSuperEnhancers
- loadPeaks
- pysam_getPeaksAt
- bedtools_getPeaksAt
- makeProfiles
- getPeaksAt
- computeMeanCov
- substractPeaks
- simpleMergePeaks
- findpeakpath
- findBestPeak
- mergeReplicatePeaks
- col_nan_scatter
- col_nan_kde_histo
- findAdditionalPeaks
- negLogLikelihood
- poissonFit
- putInBed
- pairwiseOverlap
- enrichment
- findAdditionalCobindingSignal
- fullDiffPeak
- diffPeak
- MakeSuperEnhancers
- runChromHMM
- loadMEMEmotifs
- simpleMergeMotifs
- substractPeaksTo

## highly recommended packages

*This package won't contain anything that overlap with those and might use those packages for what it is doing.*
- Bedtools
- deepTools
- MACS2
- ROSE
- MEME
- ChromHMM
