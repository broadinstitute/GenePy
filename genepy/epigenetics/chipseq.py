import os
import pandas as pd
import pysam
import numpy as np
from genepy.utils import helper as h
from genepy.utils import plot
import re
from pybedtools import BedTool
import seaborn as sns
import pyBigWig
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import minimize
from scipy.stats import zscore, fisher_exact
from scipy.special import factorial
import subprocess
import warnings
import itertools
from statsmodels.stats.multitest import multipletests

size = {"GRCh37": 2864785220, "GRCh38": 2913022398}

cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
				'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
				'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']

chroms = {'chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX',
 'chrY', '1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6',
 '7', '8', '9', 'X', 'Y'}


def bigWigFrom(bams, folder="", numthreads=8, genome='GRCh37', scaling=None, verbose=1):
	"""
	run the bigwig command line for a set of bam files in a folder

	Can apply some scaling and process files in parallel

	Args:
	----
		bams: list[str] of filepath to bam files
		folder: str folder where the bam files would be stored (will save the bigwigs there as well, else to current path ./)
		numthreadds: int number of threads to process in parallel
		genome: str genome version to use (for genome size) only GRCh37 GRCh38 (see size global variable in this file)
		scaling: list[float] a list of scaling values to apply to the bigwigs. if provided, won't scale to read counts.
		verbose: int verbose level

	"""
	if "bigwig" not in os.listdir(folder if folder else '.'):
		os.mkdir(folder + "bigwig")
	for i, bam in enumerate(bams):
		in1 = folder + bam
		out1 = folder + "bigwig/" + bam.split('/')[-1].split('.')[0] + '.bw'
		cmd = "bamCoverage --effectiveGenomeSize " + str(size[genome]) + " -p " + str(numthreads) +\
			" -b " + in1 + " -o " + out1
		if scaling is not None:
			cmd += ' --scaleFactor ' + str(scaling[i])
		if verbose == 0:
			cmd += ' 2> ' + bam + '.error.log'
		res = subprocess.run(cmd, capture_output=True, shell=True)
		if res.returncode != 0:
			raise ValueError('issue with the command: ' + str(res.stderr))


def ReadRoseSuperEnhancers(roseFolder, containsINPUT=True, name="MV411"):
	"""
	reads ROSE2's output and returns its superenhancer bedfile as a pd dataframe. 

	Can be multiple superenhancers from a set of HRK27ac bam files

	Args:
	-----
		roseFolder: str folderpath to ROSE's output
		containsINPUT: bool whether the bedfile contains INPUT signal as well
		name: str sample name from which we computed superenhancers

	Returns:
	--------
		a dataframe of the bed representation of the superenhancers
	"""
	beds = os.listdir(roseFolder)
	superenhan = pd.DataFrame()
	for i in beds:
		if i.endswith('_SuperEnhancers.table.txt'):
			superenhan = superenhan.append(pd.read_csv(roseFolder+i, sep='\t', skiprows=5)).drop(columns=['enhancerRank','isSuper'])
	data = [superenhan.columns[6]] + superenhan.columns[8:].tolist()
	if containsINPUT:
		inputd = superenhan.columns[7]
	superenhan['name'] = [i.split(name)[1].split('_')[2] for i in superenhan['REGION_ID']]
	superenhan = superenhan.rename(columns={'CONSTITUENT_SIZE':'size','REGION_ID':"id",'CHROM':'chrom','START':'start','STOP':'end',"NUM_LOCI":'num'}).replace(np.nan,0)
	superenhan['foldchange'] = superenhan[data].sum(1)/superenhan[inputd]
	superenhan = superenhan.drop(columns=data+[inputd])
	return superenhan.sort_values(by=['chrom','start','end'])


def loadPeaks(peakFile=None, peakfolder=None, isNF=True, CTFlist=[], skiprows=0):
 	"""
	loads 1 to many peak bedfile into one pandas dataframe.

	all og the peaks will be concatenated into one dataframe. this function can 
	work with jkobject|nfcore/chipseq nextflow pipelines output (isNF)

	Args:
	-----
		peakFile: str filepath to a peak bedfile
		peakfolder: str folderpath to a folder congaining beddfile
		isNF
		CTFlist
		skiprows
	"""
	if peakfolder:
		bindings = pd.DataFrame()
		for folder in os.listdir(peakfolder):
			if isNF:
				if any(tf in folder for tf in CTFlist) or not CTFlist:
					binding = pd.read_csv(peakfolder + folder + "/NA_peaks.narrowPeak", sep='\t', header=None) if\
					os.exists(peakfolder + folder + "/NA_peaks.narrowPeak") else peakfolder + folder + "/NA_peaks.broadPeak"
					binding['name'] = folder.replace('.narrowPeak', '').replace('.broadPeak','')
					bindings = bindings.append(binding)
			else:
				file = folder
				if file[-10:] in ["narrowPeak",".broadPeak"] and (any(tf in file for tf in CTFlist) or not CTFlist):
					print('reading: '+file)
					binding = pd.read_csv(peakfolder + file, sep='\t', header=None, skiprows=skiprows)
					binding['name'] = file.replace('.narrowPeak', '').replace('.broadPeak','')
					bindings = bindings.append(binding)
	elif peakFile:
		bindings = pd.read_csv(peakFile, sep='\t', header=None, skiprows=skiprows)
		bindings['name'] = peakFile.split('/')[-1].split('.')[0]
	else:
		raise ValueError("need to provide one of peakFile or peakfolder")
	if not len(bindings.columns)==6:
		bindings = bindings.drop(5, 1).drop(4, 1).rename(columns={6:4, 7:5, 8:6, 9:7,})
	bindings = bindings.rename(columns={
		0: "chrom",
		1: 'start',
		2: 'end',
		3: 'peak_number',
		4: "foldchange",
		5: "-log10pvalue",
		6: "-log10qvalue",
		7: 'relative_summit_pos'})
	bindings = bindings.sort_values(by=["chrom", "start", "end"], axis=0)
	bindings.start = bindings.start.astype('int')
	bindings.end = bindings.end.astype('int')
	if not len(bindings.columns)==6:
		bindings['relative_summit_pos'] = bindings.relative_summit_pos.astype(
			float) if 'relative_summit_pos' in bindings.columns else bindings.end - bindings.start
		bindings["-log10pvalue"] = bindings["-log10pvalue"].astype('float')
		bindings['-log10qvalue'] = bindings['-log10qvalue'].astype('float')
		loc = bindings['relative_summit_pos'].isna()
		bindings.loc[bindings[loc].index, 'relative_summit_pos'] = bindings[loc].end - bindings[loc].start
		bindings.relative_summit_pos = bindings.relative_summit_pos.astype(int)
	bindings.foldchange = bindings.foldchange.astype('float')
	bindings = bindings.reset_index(drop=True)
	return bindings


def simpleMergePeaks(peaks, window=0, totpeaknumber=0, maxp=True, mergedFold="mean"):
	"""
	simply merges bedfiles from peak callers. providing a concaneted dataframe of bed-like tables

	will recompute pvalues and foldchange from that.

	"""
	peaks = peaks.sort_values(by=['chrom', 'start','end'])
	tfs = list(set(peaks['name']))
	mergedpeaksdict = {}
	remove = []
	peaknumber = 0
	merged_bed = {
		"chrom": [peaks.iloc[0]['chrom']],
		"start": [peaks.iloc[0]['start']],
		"end": [],
		"peak_number": [peaknumber + totpeaknumber],
		"foldchange": [],
		"-log10pvalue": [],
		"-log10qvalue": [],
		"relative_summit_pos": []
	}
	foldchange = [peaks.iloc[0].get('foldchange', 1)]
	log10pvalue = [peaks.iloc[0].get('-log10pvalue', 0)]
	log10qvalue = [peaks.iloc[0].get('-log10qvalue', 0)]
	relative_summit_pos = peaks.iloc[1].get('relative_summit_pos', peaks.iloc[0]['start'])
	# computes overlap by extending a bit the window (100bp?) should be ~readsize
	prev_end = peaks.iloc[0]['end']
	prev_chrom = peaks.iloc[0]['chrom']
	tfmerged = {a: [0] for a in tfs}
	tfmerged[peaks.iloc[0]['name']][-1] = peaks.iloc[0].get('foldchange', 1)
	for i, (pos, peak) in enumerate(peaks.iloc[1:].iterrows()):
		print(str(i / len(peaks)), end="\r")
		if prev_end + window > peak['start'] and prev_chrom == peak['chrom']:
			# can be merged
			if peak.get('foldchange', 1) > max(foldchange):
				relative_summit_pos = peak.get('relative_summit_pos', peaks['start'])
			foldchange.append(peak.get('foldchange', 1))
			log10pvalue.append(peak.get('-log10pvalue', 0))
			log10qvalue.append(peak.get('-log10qvalue', 0))

		else:
			# newpeak
			for k, val in tfmerged.items():
				val.append(0)
			peaknumber += 1
			merged_bed['chrom'].append(peak['chrom'])
			merged_bed['start'].append(peak['start'])
			merged_bed['end'].append(prev_end)
			merged_bed['peak_number'].append(peaknumber + totpeaknumber)
			if mergedFold=="mean":
				merged_bed['foldchange'].append(np.mean(foldchange))
			elif mergedFold=="max":
				merged_bed['foldchange'].append(max(foldchange))
			elif mergedFold=="sum":
				merged_bed['foldchange'].append(sum(foldchange))
			else:
				raise ValueError("mergedFold needs to be one of:")
			merged_bed['-log10pvalue'].append(max(log10pvalue) if maxp else np.prod(log10pvalue))
			merged_bed['-log10qvalue'].append(max(log10qvalue) if maxp else np.prod(log10qvalue))
			merged_bed['relative_summit_pos'].append(relative_summit_pos)
			foldchange = [peak.get('foldchange', 1)]
			log10pvalue = [peak.get('-log10pvalue', 0)]
			log10qvalue = [peak.get('-log10qvalue', 0)]
			relative_summit_pos = peak.get('relative_summit_pos', peak['start'])
		prev_end = peak['end']
		prev_chrom = peak['chrom']
		tfmerged[peak['name']][-1] = peak.get('foldchange', 1)
	merged_bed['end'].append(prev_end)
	if mergedFold=="mean":
		merged_bed['foldchange'].append(np.mean(foldchange))
	elif mergedFold=="max":
		merged_bed['foldchange'].append(max(foldchange))
	elif mergedFold=="sum":
		merged_bed['foldchange'].append(sum(foldchange))
	else:
		raise ValueError("mergedFold needs to be one of:")
	merged_bed['-log10pvalue'].append(max(log10pvalue) if maxp else np.prod(log10pvalue))
	merged_bed['-log10qvalue'].append(max(log10qvalue) if maxp else np.prod(log10qvalue))
	merged_bed['relative_summit_pos'].append(relative_summit_pos)

	merged_bed = pd.DataFrame(merged_bed)
	tfmerged = pd.DataFrame(tfmerged)
	return pd.concat([merged_bed, tfmerged], axis=1, sort=False)


def findpeakpath(folder, proteiname):
	"""
	given a folder of bigwigs and a protein name, finds the right bigwig
	"""
	res = None
	for val in os.listdir(folder):
		if str(proteiname) in val:
			if res:
				raise ValueError('more than 1 bigwig file found')
			res= val
	if res:
		return res
	raise ValueError('no bigwig file found')


def findBestPeak(presence):
	"""
	given a list of -sets of peak locations for each replicate- will return the best replicate given a simple metric
	"""
	tot = []
	for ind, el in enumerate(presence):
		val = len(el)
		pres = [x for j,x in enumerate(presence) if j!=ind]
		for jnd in range(1, len(pres)+1):
			for comb in itertools.combinations(pres, jnd):
				ov = el
				for knd in range(jnd):
					ov = ov & comb[knd]
				val += len(ov)*(jnd+1)
		tot.append(val)
	return np.argsort(tot)[::-1]


def mergeReplicatePeaks(peaks, bigwigfolder, markedasbad=None, window=100,
						sampling=3000, mincov=4, doPlot=True, cov={}, minKL=8, use='max',
						MINOVERLAP=0.3, lookeverywhere=True, only='', saveloc=''):
	"""
	
	
	/!/ should only be passed peaks with at least one good replicate
	for each TFpeaksets,
	1. find the replicate that have the most peaks
	2. correlate peaks and get in highest correlation order with the replicate found in 1
	3. find overlap of both and get size of second replicate
	4. if small(er)-> use only to increase statistics
	  1. if a lot of uncalled peaks in replicate 2 at replicate 1 peaks (flag for mergebam)
	5. if similar size -> get only intersect
	  2. add to intersect, find uncalled peaks in both replicates which are called in the other
	6. repeat for all replicates
	-------------------------
	if full overlap of one of the peak replicate, only use the overlapped one to increase confidence on peak
	if >80% average non overlap,
	  print warning and percentage of overlap

	if <20% average non overlap,
	  take the overlap and increase confidence and avg logfold

	if one is <20%:
	  if other <40% average non overlap,
		take the overlap and increase confidence and avg logfold
	  else
		take

	gets the max cov at the genomic window and if above some threshold, accepts the peak.

	extend peak by X bp if no TSS
	remove TSS from peaks

	create a new data frame containing merged peak size, reassembled peak data (p value etc..) and
	a the value for presence of each TF listed in previous df
	------------------------------------

	args:
	----
	peaks: df[bed-like] all the peaks into the sameBam with a column containing the 'name'
	being the id of the sample, the 'replicate' number of this sample, the 'tf' chiped here
	bamfolder: str, foldername
	avgCov: dict(filename:int) a dict where for each bam filename is given an averageCoverage
	if use=='max':
		window:
		mincov:

	if use=='max':


	returns:
	-------
	mergedpeaks: dict{df-peakslike}
	bamtomerge: [[bam1,bam2]]

	"""
	def col_nan_scatter(x, y, **kwargs):
		df = pd.DataFrame({'x': x[:], 'y': y[:]})
		df = df[df.sum(0) != 0]
		x = df['x']
		y = df['y']
		plt.gca()
		plt.scatter(x, y)
	def col_nan_kde_histo(x, **kwargs):
		df = pd.DataFrame({'x':x[:]})
		df = df[df['x']!=0]
		x = df['x']
		plt.gca()
		sns.kdeplot(x)
	print("/!/ should only be passed peaks with at least one good replicate")
	# for a df containing a set of peaks in bed format and an additional column of different TF
	tfs = list(set(peaks['tf']))
	totpeaknumber = 0
	mergedpeaksdict = {}
	remove = []
	tomergebam = []
	ratiosofunique = {}
	h.createFoldersFor(saveloc)
	f = open(saveloc+'results.txt', 'w')
	warnings.simplefilter("ignore")
	for tf in tfs:
		if only and tf!=only:
			continue
		cpeaks = peaks[peaks.tf==tf]
		print('_____________________________________________________')
		f.write('_____________________________________________________' + '\n')
		if len(set(cpeaks['replicate'])) == 1:
			if cpeaks.name.tolist()[0] in markedasbad:
				print('the only replicate is considered bad!')
				f.write('the only replicate is considered bad!'+"\n")
				print('wrong TF: '+tf)
				f.write('wrong TF: '+tf+"\n")
				mergedpeaksdict.update({tf: cpeaks})
				remove.append(tf)
				continue
			print("we only have one replicate for " + tf + " .. pass")
			f.write("we only have one replicate for " + tf + " .. pass"+"\n")
			mergedpeaksdict.update({tf: cpeaks})
			continue
		print("merging " + tf + " peaks")
		f.write("merging " + tf + " peaks"+"\n")
		merged = simpleMergePeaks(cpeaks, window=window, maxp=False)
		merged_bed = merged[merged.columns[8:]]
		finalpeaks = merged[merged.columns[:8]]
		print('--> finish first overlaps lookup')
		f.write('--> finish first overlaps lookup'+"\n")
		# flag when  biggest is <1000 peaks
		if len(finalpeaks) < 1000:
			print('!TF has less than 1000 PEAKS!')
			f.write('!TF has less than 1000 PEAKS!'+"\n")
		# for each TF (replicates), compute number of peaks
		peakmatrix = merged_bed.values.astype(bool)

		presence = []
		for peakpres in peakmatrix.T:  # https://github.com/tctianchi/pyvenn
			presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
		# compute overlap matrix (venn?)
		if peakmatrix.shape[1] < 7 and doPlot:
			plot.venn(presence, [i+'_BAD' if i.split('-')[0] in markedasbad else i for i in merged_bed.columns], title=tf+"_before_venn", folder=saveloc)
			plt.show()
		else:
			print('too many replicates for Venn: '+str(peakmatrix.shape[1]))
			f.write('too many replicates for Venn: '+str(peakmatrix.shape[1])+"\n")
		if doPlot:
			fig = sns.pairplot(merged_bed,corner=True, diag_kind="kde", kind="reg", plot_kws={"scatter_kws":{"alpha":.05}})
			#fig = fig.map_upper(col_nan_scatter)
			#fig = fig.map_upper(col_nan_kde_histo)
			plt.suptitle("correlation of peaks in each replicate", y=1.08)
			if saveloc:
				fig.savefig(saveloc+tf+"_before_pairplot.pdf")
			plt.show()
			for i, val in enumerate(merged_bed):
				unique_inval = np.logical_and(np.delete(peakmatrix,i,axis=1).sum(1).astype(bool)==0, peakmatrix[:,i])
				sns.kdeplot(merged_bed[val][unique_inval], legend=True).set(xlim=(0,None))
			plt.title("distribution of unique peaks in each replicate")
			if saveloc:
				plt.savefig(saveloc+tf+"_before_unique_kdeplot.pdf")
			plt.show()

		bigwigs = os.listdir(bigwigfolder)

		foundgood=False
		sort = findBestPeak(presence)
		for ib,sb in enumerate(sort):
			if merged_bed.columns[sb].split('-')[0] not in markedasbad:
				foundgood=True
				break
		if not foundgood:
			print('no peaks were good enough quality')
			f.write('no peaks were good enough quality'+"\n")
			print('bad TF: '+tf)
			f.write('bad TF: '+tf+"\n")
			remove.append(tf)
			ib = 0
		# distplot
		# correlation plot


		biggest_ind = sort[ib]
		peakmatrix = peakmatrix.T
		biggest = merged_bed.columns[biggest_ind]
		print('-> main rep is: '+str(biggest))
		f.write('-> main rep is: '+str(biggest)+'\n')
		tot = peakmatrix[biggest_ind].copy().astype(int)
		# starts with highest similarity and go descending
		j = 0
		recovered = 0
		additionalpeaksinbig = np.array([])
		for i, val in enumerate(sort):
			if i==ib:
				continue
			j+=1
			# if avg non overlap > 60%, and first, and none small flag TF as unreliable.
			overlap = len(presence[biggest_ind] & presence[val]) / len(presence[biggest_ind])
			peakname = merged_bed.columns[val]
			print('- '+peakname)
			f.write('- '+peakname+'\n')
			print('  overlap: ' + str(overlap*100)+"%")
			f.write('  overlap: ' + str(overlap*100)+"%"+'\n')
			if overlap < MINOVERLAP:
				smallsupport = len(presence[biggest_ind] & presence[val]) / len(presence[val])
				print(' --> not enough overlap')
				f.write(' --> not enough overlap'+'\n')
				if smallsupport < MINOVERLAP:
					# if the secondary does not have itself the required support
					if j == 1 and merged_bed.columns[val].split('-')[0] not in markedasbad:
						print("  Wrong TF: "+tf)
						f.write("  Wrong TF: "+tf+'\n')
						remove.append(tf)
						break
					# if not first, throw the other replicate and continue
					print("  not using this replicate from the peakmatrix")
					f.write("  not using this replicate from the peakmatrix"+'\n')
					continue
			if lookeverywhere:
				tolookfor = peakmatrix[val] == 0
			else:
				tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
			# ones that we have in the Primary but not in the secondary
			additionalpeaksinsec = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(bigwigfolder, peakname), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
			if len(additionalpeaksinsec[additionalpeaksinsec>0])>0:
				sns.kdeplot(additionalpeaksinsec[additionalpeaksinsec>0],label=peakname, legend=True).set(xlim=(0,None))
				print('  min,max from newly found peaks: '+str((additionalpeaksinsec[additionalpeaksinsec>0].min(),additionalpeaksinsec[additionalpeaksinsec>0].max())))
				f.write('  min,max from newly found peaks: '+str((additionalpeaksinsec[additionalpeaksinsec>0].min(),additionalpeaksinsec[additionalpeaksinsec>0].max()))+'\n')
			# for testing purposes mainly
			finalpeaks[additionalpeaksinsec.astype(bool)].to_csv('additionalpeaksinsec_mp'+merged_bed.columns[val]+'.bed',sep='\t',index=None,header=False)
			peakmatrix[val] = np.logical_or(peakmatrix[val], additionalpeaksinsec.astype(bool))
			overlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
			if overlap < MINOVERLAP:
				newsmalloverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
				print("  we did not had enough initial overlap.")
				f.write("  we did not had enough initial overlap."+'\n')
				if newsmalloverlap < MINOVERLAP:
					if merged_bed.columns[val].split('-')[0] in markedasbad:
						print('  replicate ' + merged_bed.columns[val] + ' was too bad and had not enough overlap')
						f.write('  replicate ' + merged_bed.columns[val] + ' was too bad and had not enough overlap'+'\n')
						continue
					elif h.askif("we have two good quality peaks that don't merge well at all: "+merged_bed.columns[val] +\
							" and " +merged_bed.columns[biggest_ind]+ " can the first one be removed?:\n  \
							overlap: "+str(overlap*100)+'%\n  smalloverlap: '+str(smalloverlap*100)+'%\n  new smalloverlap: '+str(newsmalloverlap*100)+"%"):
						continue
					else:
						print("  enough from small overlaps")
						f.write("  enough from small overlaps"+'\n')
			print(' --> enough overlap')
			f.write(' --> enough overlap'+'\n')
			recovered += np.sum(additionalpeaksinsec.astype(bool))
			if merged_bed.columns[val].split('-')[0] not in markedasbad:
				tot += peakmatrix[val].astype(int)
			# ones that we have in the Primary but not in the secondary
			if not lookeverywhere or len(additionalpeaksinbig)==0:
				tolookfor = peakmatrix[biggest_ind] == 0 if lookeverywhere else np.logical_and(peakmatrix[biggest_ind]==0, peakmatrix[val])
				additionalpeaksinbig = findAdditionalPeaks(finalpeaks, tolookfor, bigwigfolder + findpeakpath(bigwigfolder, biggest), sampling=sampling, mincov=mincov, window=window, minKL=minKL, use=use)
				if len(additionalpeaksinbig[additionalpeaksinbig>0])>0:
					sns.kdeplot(additionalpeaksinbig[additionalpeaksinbig>0],label=biggest, legend=True).set(xlim=(0,None))
					print('  min,max from newly found peaks: '+str((additionalpeaksinbig[additionalpeaksinbig>0].min(),additionalpeaksinbig[additionalpeaksinbig>0].max())))
					f.write('  min,max from newly found peaks: '+str((additionalpeaksinbig[additionalpeaksinbig>0].min(),additionalpeaksinbig[additionalpeaksinbig>0].max()))+'\n')

				peakmatrix[biggest_ind] = np.logical_or(peakmatrix[biggest_ind], additionalpeaksinbig)
				tot += additionalpeaksinbig.astype(bool).astype(int)
				recovered += np.sum(additionalpeaksinbig.astype(bool))
			print('  we have recovered ' + str(recovered)+' peaks, equal to '+ str(100*recovered/np.sum(peakmatrix[biggest_ind]))+\
				'% of the peaks in main replicate')
			f.write('  we have recovered ' + str(recovered)+' peaks, equal to '+ str(100*recovered/np.sum(peakmatrix[biggest_ind]))+\
				'% of the peaks in main replicate'+'\n')
			if overlap < (MINOVERLAP+0.2)/1.2:
				# we recompute to see if the overlap changed
				newoverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[biggest_ind])
				smalloverlap = np.sum(np.logical_and(peakmatrix[val],peakmatrix[biggest_ind]))/np.sum(peakmatrix[val])
				if newoverlap < (MINOVERLAP+0.2)/1.2:
					if smalloverlap < (2+MINOVERLAP)/3:
						print("  not enough overlap to advice to merge the bams.\n  oldnew overlap: "+str(overlap*100)+'%\n  \
							new overlap: '+str(newoverlap*100)+"%")
						f.write("  not enough overlap to advice to merge the bams.\n  oldnew overlap: "+str(overlap*100)+'%\n  \
							new overlap: '+str(newoverlap*100)+"%"+'\n')
						continue
					else:
						print('  enough from small overlap to advice to merge the peaks')
						f.write('  enough from small overlap to advice to merge the peaks'+'\n')
			tomergebam.append([biggest, peakname])
			#the quality is good enough in the end we can pop from the list if it exists
			if tf in remove:
				remove.remove(tf)
		plt.title('distribution of new found peaks')
		if saveloc:
			plt.savefig(saveloc+tf+"_new_found_peaks_kdeplot.pdf")
		plt.show()
		# new distplot
		# new correlation plot
		ratiosofunique[tf] = len(np.argwhere(peakmatrix.sum(0)==1))/peakmatrix.shape[1]
		if doPlot:
			sns.pairplot(merged_bed,corner=True, diag_kind="kde", kind="reg", plot_kws={"scatter_kws":{"alpha":.05}})
			#fig = fig.map_upper(col_nan_scatter)
			#fig = fig.map_upper(col_nan_kde_histo)
			plt.suptitle("correlation and distribution of peaks after recovery", y=1.08)
			if saveloc:
				fig.savefig(saveloc+tf+"_after_pairplot.pdf")
			plt.show()
			for i, val in enumerate(merged_bed):
				unique_inval = np.logical_and(np.delete(peakmatrix,i,axis=0).sum(0).astype(bool)==0, peakmatrix[i])
				sns.kdeplot(merged_bed[val][unique_inval], legend=True).set(xlim=(0,None))
			plt.title("distribution of unique peaks in each replicate after recovery")
			if saveloc:
				plt.savefig(saveloc+tf+"_after_unique_kdeplot.pdf")
			plt.show()
		if len(peakmatrix.shape) > 1 and doPlot:
			if peakmatrix.shape[0] < 7:
				presence = []
				for peakpres in peakmatrix:  # https://github.com/tctianchi/pyvenn
					presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
				title = tf + '_recovered (TOREMOVE)' if tf in remove else tf+'_recovered'
				h.venn(presence, [i+'_BAD' if i.split('-')[0] in markedasbad else i for i in merged_bed.columns], title=title, folder=saveloc)
				plt.show()
			else:
				print('too many replicates for Venn')
				f.write('(too many replicates for Venn)'+'\n')
			finalpeaks = finalpeaks[np.logical_or(tot>1,peakmatrix[biggest_ind])]
		finalpeaks['name'] = biggest
		finalpeaks['tf'] = tf
		mergedpeaksdict.update({tf: finalpeaks})
		print(str((tf,len(finalpeaks))))
		f.write(str((tf,len(finalpeaks)))+'\n')
	mergedpeak = pd.concat([peaks for _, peaks in mergedpeaksdict.items()]).reset_index(drop=True)
	if doPlot:
		df= pd.DataFrame(data=ratiosofunique,index=['percentage of unique'])
		df['proteins'] = df.index
		fig = sns.barplot(data=df)
		plt.xticks(rotation=60,ha='right')
		plt.title("ratios of unique in replicates across experiments")
		if saveloc:
			plt.savefig(saveloc+"All_ratios_unique.pdf")
		plt.show()
	f.close()
	mergedpeak['name'] = mergedpeak.tf
	return mergedpeak, tomergebam, remove, ratiosofunique


def findAdditionalPeaks(peaks, tolookfor, filepath, sampling=1000, mincov=4,
						window=100, cov={}, minKL=8, use='max'):
	"""
	findAdditionalPeaks: for all peaks in A and/or B find in coverage file if zone has relative cov
	of more than thresh then add to peak
	if B is small and > 20% of peaks in A are found back, increase confidence and
	flag for mergeBams
	if < 20% don't flag for merge bam
	f B is big and now mean non overlap < 40%, take union and flag for mergeBam else, throw B.

	Args:
	-----
		peaks
		tolookfor
		filepath
		sampling
		mincov
		window
		cov
		minKL
		use
	returns:
	-------
		np.array(bool) for each peaks in peakset, returns a binary
	"""
	# def poisson(k, lamb, scale): return scale * (lamb**k / factorial(k)) * np.exp(-lamb)

	def KLpoisson(lamb1, lamb2): return lamb1 * np.log(lamb1 / lamb2) + lamb2 - lamb1
	def poisson(k, lamb): return (lamb**k/factorial(k)) * np.exp(-lamb)
	def negLogLikelihood(params, data): return - np.sum(np.log(poisson(data, params[0])))
	def poissonFit(data): return float(minimize(negLogLikelihood,x0=np.ones(1),args=(data,),method='Powell').x)
	bw = pyBigWig.open(filepath)
	res = np.zeros(len(peaks))
	prevchrom = ''
	lamb = {}
	cov = {}
	#ignore by message
	warnings.filterwarnings("ignore", message="encountered in")
	for i, has in enumerate(tolookfor):
		if has:
			val = peaks.iloc[i]
			if val.chrom not in chroms:
				continue
			if val.chrom != prevchrom:
				if val.chrom not in cov:
					cov[val.chrom] = bw.stats(str(val.chrom))[0]
					prevchrom = val.chrom
					if use == 'poisson':
						#TODO: compute on INPUT file instead
						samples = np.zeros(window * sampling)
						sam = np.random.rand(sampling)
						sam = sam * (bw.chroms(str(val.chrom))-window)
						for j, sample in enumerate(sam.astype(int)):
							samples[j*window:(j + 1)*window] = np.nan_to_num(bw.values(str(val.chrom), sample, sample + window), 0)
						scale = np.unique(samples)[1]
						samples = (samples/scale).astype(int)
						lamb[val.chrom] = (poissonFit(samples),scale)

			start = max([val.start - window, 0])
			end = min(val.end + window, bw.chroms(str(val.chrom)))
			zone = np.nan_to_num(bw.values(str(val.chrom), start, end), 0)
			if use == 'max':
				if max(zone) / cov[val.chrom] > mincov*1.5 or sum(zone) / (cov[val.chrom] * (end - start)) > mincov:
					res[i] = max(zone) / cov[val.chrom]
			elif use == 'poisson':
				#TODO: compute -log10pvalue
				la = poissonFit((zone/lamb[val.chrom][1]).astype(int))
				kl=KLpoisson(la, lamb[val.chrom][0])
				if kl > minKL:
					res[i] = max(zone) / cov[val.chrom] #foldchange from macs2

	return res


def putInBed(conscensus, value, window=10, mergetype='mean'):
	""" given ordered dataframes

	conscensus df[start,end,chrom]
	value df[start, end, chrom,foldchange]

	type: oneof: mean,first,last,
	"""
	conscensus = conscensus.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	value = value.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	locinvalue=0
	loc=0
	tot=0
	num = []
	res = np.zeros(len(conscensus))
	not_end=True
	def add(res,num,not_end,loc):
		if len(num)>0:
			if mergetype=='mean':
				res[loc] = np.mean(num)
			elif mergetype=='first':
				res[loc]=np.sum(num)
			elif mergetype=='first':
				res[loc]=num[0]
			elif mergetype=='last':
				res[loc]=num[-1]
			else:
				raise ValueError('must be one of')
			num= []
		loc+=1
		if loc == len(conscensus):
			not_end=False
		return res, num, not_end,loc
	while not_end:
		print(loc/len(conscensus),end="\r")
		a = conscensus.iloc[loc]
		b = value.iloc[locinvalue]
		if b.chrom < a.chrom:
			locinvalue+=1
			if locinvalue == len(value):
				not_end=False
		elif b.chrom > a.chrom:
			loc+=1
			if loc == len(conscensus):
				not_end=False
		elif b.start<a.start:
			if b.end+window>a.start:
				tot+=1
				num.append(b.foldchange)
				if b.end>a.end+window:
					res,num,not_end,loc = add(res,num,not_end,loc)
					continue
			locinvalue+=1
			if locinvalue == len(value):
				not_end=False
		elif b.start<a.end+window:
			tot+=1
			num.append(b.foldchange)
			if b.end>a.end+window:
				res,num,not_end,loc = add(res,num,not_end,loc)
				continue
			locinvalue+=1
			if locinvalue == len(value):
				not_end=False
		else:
			res,num,not_end,loc = add(res,num,not_end,loc)
	print(str(tot)+' were merged into conscensus')
	return res


def pairwiseOverlap(bedfile, norm=True, bedcol=8, correct=True, docorrelation=True):
	"""
	considering a befile representing a conscensus set of peaks
	with each columns after the 7th one representing the signal of a given ChIP experiment
	over this conscensus

	overlap of j in i
	overlap of row values in col values
	"""
	if correct:
		print("we will be correcting for fully similar lines/ columns by removing 1 on their last value")
	dat = bedfile[bedfile.columns[bedcol:]].values
	prob = dat.astype(bool).sum(0)/len(dat)
	correlation = np.ones((dat.shape[1],dat.shape[1]))
	overlap = np.ones((dat.shape[1],dat.shape[1]))
	for i, col in enumerate(dat.T):
		#pdb.set_trace()
		overlapping = np.delete(dat,i,axis=1)[col!=0]
		col = col[col!=0]
		add=0
		for j, val in enumerate(overlapping.T):
			if j==i:
				add=1
			if docorrelation:
				if norm and not np.isnan(zscore(col[val!=0]).sum()) and not np.isnan(zscore(val[val!=0]).sum()) or not correct:
					correlation[i,j+add] = np.corrcoef(zscore(val[val!=0]),zscore(col)[val!=0])[0,1]
				else:
					tmp = np.corrcoef(val[val != 0], col[val != 0])[0, 1]
					if np.isnan(tmp) and correct:
						if len(col[val!=0]) == 0 or len(val[val!=0]) == 0:
						# one has no overlap
							correlation[i,j+add] =  0
						else:
							# one contains only the same value everywhere
							col[-1]-=max(0.01,abs(np.mean(col)))
							val[-1]-=max(0.01,abs(np.mean(val)))
							correlation[i,j+add] = np.corrcoef(val,col)[0,1]
			overlap[i,j+add]=len(val[val!=0])/len(col)
	if docorrelation:
		correlation = pd.DataFrame(data=correlation, index=bedfile.columns[bedcol:], columns=bedfile.columns[bedcol:])
		correlation[correlation.isna()] = 0
	overlap = pd.DataFrame(data=overlap, index=bedfile.columns[bedcol:], columns=bedfile.columns[bedcol:]).T
	return overlap, correlation if docorrelation else None


def enrichment(bedfile, bedcol=8, groups=None, okpval=10**-3):
	"""
	considering a befile representing a conscensus set of peaks
	with each columns after the 7th one representing the signal of a given ChIP experiment
	over this conscensus

	enrichment of j in i
	enrichment of row values in col values
	"""
	dat = bedfile[bedfile.columns[bedcol:]].values
	# pdb.set_trace()
	prob = dat.astype(bool).sum(0)/len(dat)
	enrichment = np.zeros((dat.shape[1] if groups is None else len(set(groups)), dat.shape[1]))
	pvals = np.zeros(
		(dat.shape[1] if groups is None else len(set(groups)), dat.shape[1]))
	if groups is not None:
		for i in set(groups):
			overlapping = dat[groups==i]
			for j,val in enumerate(overlapping.T):
				# enrichment of j in i
				e, p = fisher_exact([
					[len(val[val != 0]), len(val[val == 0])],
					[prob[j]*len(dat), (1-prob[j])*len(dat)]])
				enrichment[i, j] = np.log2(e)
				pvals[i, j] = p
	else:
		for i, col in enumerate(dat.T):
			overlapping = np.delete(dat,i,axis=1)[col!=0]
			col = col[col!=0]
			add=0
			for j, val in enumerate(overlapping.T):
				if j==i:
					add=1
					enrichment[i,i]=0
				e, p = fisher_exact([[len(val[val != 0]), len(val[val == 0])], [
					prob[j+add]*len(dat), (1-prob[j+add])*len(dat)]])
				enrichment[i, j+add] = np.log2(e)
				pvals[i, j+add] = p
		enrichment[i,i]=0
	enrichment = pd.DataFrame(data=enrichment, index=bedfile.columns[bedcol:] if groups is None else set(groups), columns=bedfile.columns[bedcol:]).T
	enrichment[enrichment==-np.inf] = -1000
	enrichment[enrichment.isna()] = 0
	enrichment[enrichment == np.inf] = 1000
	pvals = np.reshape(multipletests(pvals.ravel(),
									 0.1, method="bonferroni")[1], pvals.shape)
	pvals = pd.DataFrame(
		data=pvals, index=bedfile.columns[bedcol:]  if groups is None else set(groups), columns=bedfile.columns[bedcol:]).T
	enrichment[pvals>okpval] = 0
	return enrichment, pvals


def findAdditionalCobindingSignal(conscensus, known=None, bigwigs=[], window=100):
	"""
	somewhat similar concept to computePeaksAt

	# get pysam data
	# ask for counts only at specific locus based on peak windows from mergedpeakset
	# append to an array
	# return array, normalized
	"""
	if known:
		print('getting '+ str(len(peaks.tf))+' peaks. Using the peaks values directly if \
				available and using the bigwigs otherwise.')
		res = known.values.astype(float)
	elif len(bigwigs)>0:
		print('getting '+str(len(bigwigs))+' bigwigs, no peaks passed. Will compute the cobinding values\
			across the conscensus for each bigwigs.')
		res = np.zeros((len(bigwigs), len(conscensus)), dtype=float)
	else:
		raise ValueError('you need to pass a list of path to bigwigs for each/some samples')
	for i, bw in enumerate(bigwigs):
		if known:
			found = False
			for j, val in enumerate(known.tf):
				if val in bw:
					if found:
						raise ValueError('found two or more matching tf for bigwig: '+str(bw))
					found = True
					i=j
					break
				if not found:
					print('no tf found in known for tf: '+bw)
					raise ValueError('you need to have an amount of known columns equal to your bigwigs')
		print('doing file ' + str(bw))
		bw = pyBigWig.open(bw)
		for k, val in conscensus.iterrows():
			if known:
				if res[i][k]!=0:
					continue
			start = max([val.start - window, 0])
			end = min(val.end + window, bw.chroms(str(val.chrom)))
			res[i][k] = bw.stats(str(val.chrom), start, end)[0]
	res = np.nan_to_num(res, 0)
	return conscensus.join(pd.Dataframe(data=(res.T / res.max(1)).T, columns=bigwigs if not known else known.columns))



def annotatePeaks():
	"""
	get H3k27AC peaks and compute Rose
	for each TF peak
	assign super enhancer if within it and create a super enhancer TFgroup
	for each peaks
	apply similar to rose and merge TF together. (create new TFgroup)
	for each TF groups say
	  if within super enhancer
	  if within high h3k27ac marks
	  its nearest gene (distance to it, cis/trans)


	FIND A WAY TO FILTER TO ONLY PLACES WITH H3K27AC marks??

	TAD points where most contacts on one side happened on this side.
	specific distance. zone of most concentration of contacts for a peak region.
	"""
	# if way == 'Closest':

	# elif way == 'ClosestExpressed':

	# elif way == 'ActivityByContact':

	# else:
	#     raise ValueError('needs to be oneof Closest ClosestExpressed ActivityByContact')


# os.system()


def getCoLocalization():
	"""
	for each annotations (super enhancer & TFgroups)
	for each TF, find highest peak/meanCov. if above thresh add to localization
	"""


def refineGroupsWithHiC():
	"""
	given HiC data, for each loops (should be less than X Mb distance. create a Xkb zone around both sides
	and find + merge each TF/TFgroup at this location)
	"""


def fullDiffPeak(bam1, bam2, control1, size=None, control2=None, scaling=None, directory='diffData/',
				 res_directory="diffPeaks/", isTF=False, compute_size=True, pairedend=True):
	"""
	will use macs2 to call differential peak binding from two bam files and their control

	one can also provide some spike in scaling information

	Args:
	-----
	bam1
	bam2()
	control1
	control2
	scaling
	"""
	print("doing diff from " + bam1 + " and " + bam2)
	if scaling is not None:
		if max(scaling) > 1:
			raise ValueError("scalings need to be between 0-1")
	name1 = bam1.split('/')[-1].split('.')[0]
	name2 = bam2.split('/')[-1].split('.')[0]
	if size is None:
		if isTF:
			size = 147
		else:
			size = 200
	if compute_size:
		print('computing the fragment avg size')
		cmd = "macs2 predictd -i " + bam1
		ret = subprocess.run(cmd, capture_output=True, shell=True)
		size = re.findall("# predicted fragment length is (\d+)", str(ret.stderr))[0]
		print(size)
	else:
		print('using default|given size')
	pairedend = "BAMPE" if pairedend else "BAM"
	if control2 is None:
		control2 = control1
	cmd1 = "macs2 callpeak -B -t " + bam1 + " -c " + control1 + " --nomodel --extsize " + str(size) + " -n " + name1 + " --outdir " + directory + " -f " + pairedend
	cmd2 = "macs2 callpeak -B -t " + bam2 + " -c " + control2 + " --nomodel --extsize " + str(size) + " -n " + name2 + " --outdir " + directory + " -f " + pairedend
	print('computing the scaling values')
	ret = subprocess.run(cmd1, capture_output=True, shell=True)
	print(ret.stderr)
	scaling1a = int(re.findall(" after filtering in treatment: (\d+)", str(ret.stderr))[0])
	scaling1b = int(re.findall(" after filtering in control: (\d+)", str(ret.stderr))[0])
	scaling1 = scaling1a if scaling1a <= scaling1b else scaling1b
	ret = subprocess.run(cmd2, capture_output=True, shell=True)
	print(ret.stderr)
	scaling2a = int(re.findall(" after filtering in treatment: (\d+)", str(ret.stderr))[0])
	scaling2b = int(re.findall(" after filtering in control: (\d+)", str(ret.stderr))[0])
	scaling2 = scaling2a if scaling2a <= scaling2b else scaling2b
	if scaling is not None:
		scaling1 = int(scaling1/scaling[0])
		scaling2 = int(scaling2/scaling[1])
	print(scaling1, scaling2)
	return diffPeak(directory+name1+"_treat_pileup.bdg", directory+name2+"_treat_pileup.bdg",
		directory+name1+"_control_lambda.bdg", directory+name2+"_control_lambda.bdg",
		res_directory, scaling1, scaling2, size)


def diffPeak(name1, name2, control1, control2, res_directory, scaling1, scaling2, size):
	"""
	calls MACS2 bdgdiff given the parameters
	"""
	print("doing differential peak binding")
	cmd = "macs2 bdgdiff --t1 " + name1 + " --c1 "
	cmd += control1+" --t2 " + name2 +" --c2 " + control2
	cmd += " --d1 " + str(scaling1) + " --d2 " + str(scaling2) + " -g 60 "
	cmd += "-l " + str(size) + " --o-prefix " + name1.split('/')[-1].split('.')[0] + "_vs_"
	cmd += name2.split('/')[-1].split('.')[0] + " --outdir " + res_directory
	res = subprocess.run(cmd, capture_output=True, shell=True)
	return res


def AssignToClosestExpressed(bed,countFile,genelocFile):
	print("the bed file and genelocFile should use the same assembly")
	genelocFile = pd.read_csv(genelocFile,sep="\t", compression="", columns=['chrom','start','end',"name"])
	#for val in bed.iterrows():


def MakeSuperEnhancers(MACS2bed, bamFile, outdir, baiFile=None, rosePath=".",
	stitching_distance=None, TSS_EXCLUSION_ZONE_SIZE="2500", assembly="hg38",controlBam=None,controlBai=None):
	"""
	Calls super enhancer from H3K27ac with the ROSE algorithm

	Args:
	----
		MACS2GFF
		bamFile
		outdir
		baiFile
		rosePath
		stitching_distance
		TSS_EXCLUSION_ZONE_SIZE
		assembly
		controlBam
		controlBai

	Returns:
	--------
		a bed-like dataframe with the superenhancers

	outdir has to be an absolute path or a ~/path
	"""
	print("we are going to move your input files to "+rosePath)
	cmd = 'mv '+MACS2bed+' '+rosePath+MACS2bed.split('/')[-1]+'.bed'
	cmd += ' && mv '+bamFile+' '+rosePath
	baiFile = baiFile if baiFile else bamFile[:-1]+'i'
	cmd += ' && mv '+baiFile +' '+rosePath
	res = subprocess.run(cmd, capture_output=True, shell=True)
	if res.returncode != 0:
		raise SystemError('failed moving files: '+str(res))
	cmd = "cd "+rosePath+" && python ROSE_main.py -g "+assembly+" -i "+MACS2bed.split('/')[-1]+'.bed' + " -r " + bamFile.split('/')[-1] + " -o " + outdir
	if TSS_EXCLUSION_ZONE_SIZE:
		cmd+=" -t "+TSS_EXCLUSION_ZONE_SIZE
	if controlBam:
		os.system('mv '+controlBam+' '+rosePath)
		controlBai = controlBai if controlBai else controlBam[:-1]+'i'
		os.system('mv '+controlBai+' '+rosePath)
		cmd+=" -c "+controlBam.split('/')[-1]
	if stitching_distance:
		cmd+=" -s "+ stitching_distance

	res = subprocess.run(cmd, capture_output=True, shell=True)
	fail = False
	if res.returncode != 0:
		v = 'ROSE failed:' +str(res)
		fail = True
	print('finished.. moving them back to their original folder')
	cmd = 'mv '+rosePath+MACS2bed.split('/')[-1]+'.bed ' + MACS2bed
	cmd += ' && mv '+rosePath+bamFile.split('/')[-1]+' '+bamFile
	cmd+= ' && mv '+rosePath+baiFile.split('/')[-1]+' '+baiFile
	if controlBam:
		cmd += ' && mv '+rosePath+controlBam.split('/')[-1]+' '+controlBam
		cmd += ' && mv '+rosePath+controlBai.split('/')[-1]+' '+controlBai
	res = subprocess.run(cmd, capture_output=True, shell=True)
	if res.returncode != 0:
		raise SystemError('failed moving files: '+str(res))
	if fail:
		raise SystemError(v)
	print('worked')
	print(res)
	return ReadRoseSuperEnhancers(outdir ,bool(controlBam))



def runChromHMM(outdir, data, numstates=15, datatype='bed', folderPath=".", chromHMMFolderpath="~/ChromHMM/", assembly="hg38",control_bam_dir=None):
	"""
	runs chromHMM algorithm

	Args:
	-----
		outdir str: an existing dir where the results should be saved
		data: a df[cellname,markname,markbed|bam|bigwig, ?controlbed|bam|bigwig]
		folderpath
		numstates
		datatype
		folderPath
		chromHMMFolderpath
		assembly
		control_bam_dir

	Returns:
	-------
		A dict of bed like dataframes containing the regions of the different states
	"""
	print("you need to have ChromHMM")
	chromHMM = "java -mx8000M -jar "+chromHMMFolderpath+"ChromHMM.jar "
	h.createFoldersFor(outdir+'binarized/')
	data.to_csv(outdir+"input_data.tsv", sep='\t',index=None,header=None)
	cmd = chromHMM
	if datatype=="bed":
		cmd+="BinarizeBed "
	elif datatype=="bigwig":
		cmd+="BinarizeSignal "
	elif datatype=="bam":
		cmd+="BinarizeBam "
	else:
		raise ValueError('you need to provide one of bam, bigwig, bed')
	cmd+= chromHMMFolderpath+"CHROMSIZES/"+assembly+".txt "+ folderPath+" "+outdir+"input_data.tsv "+outdir+"binarized"
	if control_bam_dir:
		cmd+=" -c "+control_bam_dir
	res1 = subprocess.run(cmd, capture_output=True, shell=True)
	print(res1)
	if res1.returncode!=0:
		raise ValueError(str(res1.stderr))
	cmd = chromHMM + "LearnModel -printposterior -noautoopen "
	if len(data)<10:
		cmd += '-init load -m '+chromHMMFolderpath+'model_15_coreMarks.txt '
	cmd += outdir+"binarized "+outdir+" "+str(numstates)+" "+assembly
	res2 = subprocess.run(cmd, capture_output=True, shell=True)
	print(res2)
	if res2.returncode!=0:
		raise ValueError(res2.stderr)
	ret = {}
	for v in set(data[0]):
		ret[v] = pd.read_csv(outdir+v+'_'+str(numstates)+'_dense.bed', sep='\t', header=None,
			skiprows=1).drop(columns=[4,5,6,7]).rename(columns=
			{0:'chrom',1:'start',2:'end',3:'state',8:"color"})
	return ret

def loadMEMEmotifs(file, tfsubset=[],motifspecies='HUMAN'):
	if file.endswith('.gff'):
		print('converting to bed, you need to have "gfftobed" installed')
		cmd = 'gff2bed < '+file+' > '+file+'.bed'
		file = file+'.bed'
		res = subprocess.run(cmd,capture_output=True, shell=True)
		if res.returncode != 0:
			raise ValueError('issue with the command: ' + str(res.stderr))
		else:
			print(res.stdout.decode("utf-8"))
	## What are the motifs of our CRC members in ATACseq but not in our matrix
	merged_motif = pd.read_csv(file, sep='\t',skiprows=0,index_col=None, names=['pos',"fimo",
		"nucleotide_motif","relStart","relEnd","pval","strand",".","data"])
	merged_motif['tf']=[i[5:].split("_"+motifspecies)[0] for i in merged_motif.data]
	if tfsubset:
		merged_motif = merged_motif[merged_motif.tf.isin(tfsubset)]
	merged_motif['chrom'] = [i.split(':')[0][3:] for i in merged_motif.index]
	merged_motif['start'] = [i.split(':')[1].split('-')[0] for i in merged_motif.index]
	merged_motif['end'] = [i.split(':')[1].split('-')[1] for i in merged_motif.index]
	merged_motif = merged_motif.reset_index(drop=True)
	merged_motif = merged_motif.rename(columns={'pos':'relStart','fimo':'relEnd','nucleotide_motif':'pos',
		'relStart':'pval','relEnd':'strand','pval':'fimo','strand':'motif'})
	merged_motif['motif'] = [i.split('sequence=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif['p_val'] = [i.split('pvalue=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif['q_val'] = [i.split('qvalue=')[1].split(';')[0] for i in merged_motif.data]
	merged_motif = merged_motif.drop(columns=['pos','.','fimo','data'])
	merged_motif = merged_motif[merged_motif.columns[[6,7,8,0,1,3,2,9,10,5,4]]]
	merged_motif = merged_motif.sort_values(by=['chrom','start','end']).reset_index(drop=True)
	return merged_motif


def simpleMergeMotifs(motifs, window=0):
	if type(motifs) is list:
		motifs = pd.concat(motifs)
	motifs = motifs.sort_values(by=['chrom', 'start'])
	toremove = []
	issues = []
	prevmotif = motifs.iloc[0]
	for i, (pos, motif) in enumerate(motifs.iloc[1:].iterrows()):
		print(str(i / len(motifs)), end="\r")
		if prevmotif['end'] + window > motif['start'] and prevmotif['chrom'] == motif['chrom']:
			# can be merged
			if motif['tf']!= prevmotif['tf'] or motif['motif'] != prevmotif['motif']:
				print('found different motifs overlapping')
				issues.extend([motif,prevmotif])
			else:
				toremove.append(pos)
		prevmotif = motif
	motifs = motifs.drop(index=toremove).reset_index(drop=True)
	issues = pd.concat(issues)
	return motifs, issues


def substractPeaksTo(peaks,loci, bp=50):
	"""
	removes all peaks that are not within a bp distance to a set of loci

	Args:
	----
		peaks: a bed file df with a chrom,start, end column at least
		loci: a df witth a chrom & loci column
		bp: the max allowed distance to the loci

	Returns:
	-------
		all the peaks that are within this distance
	"""
	i=0
	j=0
	keep=[]
	bp=50
	while j<len(peaks) and i<len(loci):
		h.showcount(j,len(peaks))
		if peaks.loc[j].chrom > loci.loc[i].chrom:
			i+=1
			continue
		if peaks.loc[j].chrom < loci.loc[i].chrom:
			j+=1
			continue
		if peaks.loc[j].start - bp > loci.loc[i].loci:
			i+=1
			continue
		if peaks.loc[j].end + bp< loci.loc[i].loci:
			j+=1
			continue
		if peaks.loc[j].end + bp >= loci.loc[i].loci and peaks.loc[j].start - bp <= loci.loc[i].loci:
			keep.append(j)
			j+=1
	return peaks.loc[set(keep)]

	
