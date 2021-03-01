import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from genepy.epigenetics import chipseq as chip

def plotAverageOfSamples(samples, folder="", showAll=False, maxv=None, minv=None):
  res = [] 
  plt.figure()
  plt.ylim(minv,maxv)
  for sample in samples:
    data = pd.read_csv(sample, sep='\t', skiprows=1, header=None, names=['chr', 'start', 'end', 'name', "foldchange","."]+list(range(600)))
    r = data[list(range(600))].mean().tolist()
    res.append(r)
    if showAll:
      sns.lineplot(data=np.array(r), color="#BFBFFF")
  sns.lineplot(data=np.array(res).mean(0))
  if folder:
    plt.savefig(folder+"_averageofsamples.pdf", color="#1F1FFF")
  return res


def pysam_getPeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):

	# get pysam data
	# ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
	# for each counts, do a rolling average (or a convolving of the data) with numpy
	# append to an array
	# return array, normalized
	loaded = {}
	res = {i: np.zeros((len(peaks), window * 2)) for i in bams}
	peaks = peaks.sort_values(by="foldchange", ascending=False).iloc[:numpeaks]
	peaks.chrom = peaks.chrom.astype(str)
	for val in bams:
	loaded.update({val: pysam.AlignmentFile(
            folder + val, 'rb', threads=numthreads)})
	for k, bam in loaded.items():
		for num, (i, val) in enumerate(peaks.iterrows()):
			print(int(num / len(peaks)), end='\r')
			center = int((val['start'] + val['end']) / 2)
			for pileupcolumn in bam.pileup(val['chrom'], start=center - window,
                                  stop=center + window, truncate=True):
				res[k][num][pileupcolumn.pos - (center - window)] = pileupcolumn.n
	fig, ax = plt.subplots(1, len(res))
	for i, (k, val) in enumerate(res.items()):
		sns.heatmap(val, ax=ax[i])
		ax[i].set_title(k.split('.')[0])
	fig.show()
	return res, fig


def bedtools_getPeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):
	"""
	get pysam data
	ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
	for each counts, do a rolling average (or a convolving of the data) with numpy
	append to an array
	return array, normalized
	"""
	loaded = {}
	center = [int((val['start'] + val['end']) / 2) for k, val in peaks.iterrows()]
	peaks['start'] = [c - window for c in center]
	peaks['end'] = [c + window - 1 for c in center]
	peaks[peaks.columns[:3]].sort_values(by=['chrom', 'start']).to_csv(
            'temp/peaks.bed', sep='\t', index=False, header=False)
	bedpeaks = BedTool('temp/peaks.bed')

	fig, ax = plt.subplots(1, len(bams))
	peakset = peaks["foldchange"].values.argsort()[::-1][:numpeaks]
	for i, val in enumerate(bams):
		coverage = BedTool(folder + val).intersect(bedpeaks).genome_coverage(bga=True, split=True)\
			.intersect(bedpeaks).to_dataframe(names=['chrom', 'start', 'end', 'coverage'])
		cov = np.zeros((len(peaks), window * 2), dtype=int)
		j = 0
		pdb.set_trace()
		for i, (k, val) in enumerate(peaks.iterrows()):
			print(i / len(peaks), end='\r')
			while coverage.iloc[j].start > val.start:
				j -= 1
			while coverage.iloc[j].start < val.end:
				cov[i][coverage.iloc[j].start - val.start:coverage.iloc[j].end - val.start] =\
					coverage.iloc[j].coverage
				j += 1
		sns.heatmap(coverage, ax=ax[i])
		ax[i].set_title(val.split('.')[0])
	fig.show()
	return None, fig


def makeProfiles(matx=[], folder='', matnames=[], title='',
                 name='temp/peaksat.pdf', refpoint="TSS", scale=None,
                 sort=False, withDeeptools=True, cluster=1, vmax=None, vmin=None, overlap=False,
                 legendLoc=None):
	if withDeeptools:
		if not (len(matnames) == 2 and len(matx) == 2):
			raise ValueError('you need two mat.gz files and two names')
		h.createFoldersFor(name)
		cmd = 'computeMatrixOperations relabel -m '
		cmd += matx[0] + ' -o '+matx[0]+' --groupLabels '+matnames[0]
		cmd += ' && computeMatrixOperations relabel -m '
		cmd += matx[1] + ' -o '+matx[1]+' --groupLabels '+matnames[1]
		cmd += ' && computeMatrixOperations rbind -m '
		cmd += matx[0] + ' ' + matx[1] + " -o " + \
                    '.'.join(name.split('.')[:-1]) + ".gz"
		cmd += ' && plotProfile'
		cmd += " --matrixFile " + '.'.join(name.split('.')[:-1]) + ".gz"
		cmd += " --outFileName " + name
		cmd += " --refPointLabel " + refpoint
		if vmax is not None:
			cmd += " -max "+str(vmax)
		if vmin is not None:
			cmd += " -min "+str(vmin)
		if cluster > 1:
			cmd += " --perGroup --kmeans "+str(cluster)
		if legendLoc:
			cmd += " --legendLocation "+legendLoc
		if title:
			cmd += " --plotTitle " + title
		data = subprocess.run(cmd, shell=True, capture_output=True)
		print(data)


def getPeaksAt(peaks, bigwigs, folder='', bigwignames=[], peaknames=[], window=1000, title='', numpeaks=4000, numthreads=8,
               width=5, length=10, torecompute=False, name='temp/peaksat.pdf', refpoint="TSS", scale=None,
               sort=False, withDeeptools=True, onlyProfile=False, cluster=1, vmax=None, vmin=None, overlap=False,
               legendLoc=None):
	"""
	get pysam data
	ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
	for each counts, do a rolling average (or a convolving of the data) with numpy
	append to an array
	return array, normalized
	"""
	if withDeeptools:
		if isinstance(peaks, pd.DataFrame):
			peaks = 'peaks.bed '
			peaks.to_csv('peaks.bed', sep='\t', index=False, header=False)
		elif type(peaks) == list:
			pe = ''
			i = 0
			for n, p in enumerate(peaks):
				if 20 < int(os.popen('wc -l ' + p).read().split(' ')[0]):
					pe += p + ' '
				elif len(peaknames) > 0:
					peaknames.pop(n-i)
					i += 1
			peaks = pe
		elif type(peaks) == str:
			peaks += ' '
		else:
			raise ValueError(' we dont know this filetype')
		if type(bigwigs) is list:
			pe = ''
			for val in bigwigs:
				pe += folder + val + ' '
			bigwigs = pe
		else:
			bigwigs = folder + bigwigs + ' '
		h.createFoldersFor(name)
		cmd = ''
		if not os.path.exists('.'.join(name.split('.')[:-1]) + ".gz") or torecompute:
			cmd += "computeMatrix reference-point -S "
			cmd += bigwigs
			cmd += " --referencePoint "+refpoint
			cmd += " --regionsFileName " + peaks
			cmd += " --missingDataAsZero"
			cmd += " --outFileName " + '.'.join(name.split('.')[:-1]) + ".gz"
			cmd += " --upstream " + str(window) + " --downstream " + str(window)
			cmd += " --numberOfProcessors " + str(numthreads) + ' && '
		cmd += "plotHeatmap" if not onlyProfile else 'plotProfile'
		if type(name) is list:
			if not onlyProfile:
				raise ValueError('needs to be set to True, can\'t average heatmaps')
			cmd += " --matrixFile " + '.gz '.join(name) + ".gz"
			if average:
				cmd += "--averageType mean"
		else:
			cmd += " --matrixFile " + '.'.join(name.split('.')[:-1]) + ".gz"
		cmd += " --outFileName " + name
		cmd += " --refPointLabel " + refpoint
		if vmax is not None:
			cmd += " -max "+str(vmax)
		if vmin is not None:
			cmd += " -min "+str(vmin)
		if cluster > 1:
			cmd += " --perGroup --kmeans "+str(cluster)
		if overlap:
			if onlyProfile:
				cmd += " --plotType overlapped_lines"
			else:
				raise ValueError("overlap only works when onlyProfile is set")
		if legendLoc:
			cmd += " --legendLocation "+legendLoc

		if len(peaknames) > 0:
			pe = ''
			for i in peaknames:
				pe += ' ' + i
			cmd += " --regionsLabel" + pe
		if type(bigwigs) is list:
			if len(bigwignames) > 0:
				pe = ''
				for i in bigwignames:
					pe += ' "' + i + '"'
				cmd += " --samplesLabel" + pe
		if title:
			cmd += " --plotTitle '"+title+"'"
		data = subprocess.run(cmd, shell=True, capture_output=True)
		print(data)
	else:
		if 'relative_summit_pos' in peaks.columns:
		center = [int((val['start'] + val['relative_summit_pos']))
                    for k, val in peaks.iterrows()]
		else:
		center = [int((val['start'] + val['end']) / 2)
                    for k, val in peaks.iterrows()]
		pd.set_option('mode.chained_assignment', None)
		peaks['start'] = [c - window for c in center]
		peaks['end'] = [c + window for c in center]
		fig, ax = plt.subplots(1, len(bigwigs), figsize=[
		                       width, length], title=title if title else 'Chip Heatmap')
		if sort:
			peaks = peaks.sort_values(by=["foldchange"], ascending=False)
		if numpeaks > len(peaks):
			numpeaks = len(peaks) - 1
		cov = {}
		maxs = []
		for num, bigwig in enumerate(bigwigs):
			bw = pyBigWig.open(folder + bigwig)
			co = np.zeros((numpeaks, window * 2), dtype=int)
			scale = scale[bigwig] if scale is dict else 1
			for i, (k, val) in enumerate(peaks.iloc[:numpeaks].iterrows()):
				try:
					co[i] = np.nan_to_num(bw.values(str(val.chrom), val.start, val.end), 0)
				except RuntimeError as e:
					print(str(val.chrom), val.start, val.end)
					pass
			cov[bigwig] = co
			maxs.append(co.max())
		for num, bigwig in enumerate(bigwigs):
		sns.heatmap(cov[bigwig] * scale, ax=ax[num], vmax=max(maxs), yticklabels=[], cmap=cmaps[num],
                    cbar=True)
		ax[num].set_title(bigwig.split('.')[0])
		fig.subplots_adjust(wspace=0.1)
		fig.show()
		fig.savefig(name)
		return cov, fig


def andrew(groups, merged, annot, enr=None, pvals=None, cols=8, precise=True, title = "sorted clustermap of cobindings clustered", folder="", rangeval=4, okpval=10**-3, size=(20,15),vmax=3, vmin=0):
    if enr is None or pvals is None:
        enr, pvals = chip.enrichment(merged, groups=groups)
    rand = np.random.choice(merged.index,5000)
    subgroups = groups[rand]
    sorting = np.argsort(subgroups)
    redblue = cm.get_cmap('RdBu_r',256)
    subenr = enr.iloc[annot-cols:]
    subenr[subenr>rangeval]=rangeval
    subenr[subenr<-rangeval]=-rangeval
    subenr = subenr/rangeval
    data = []
    #colors = []
    impv = pvals.values
    for i in subgroups[sorting]:
        #colors.append(viridis(i))
        a = redblue((128+(subenr[i]*128)).astype(int)).tolist()
        for j in range(len(a)):
            a[j] = [1.,1.,1.,1.] if impv[j,i] > okpval else a[j]
        data.append(a)
    data = pd.DataFrame(data=data,columns=list(subenr.index),index= rand[sorting])
    #data["clusters"]  = colors
    
    a = np.log2(1.01+merged[merged.columns[cols:annot]].iloc[rand].iloc[sorting].T)
    if not precise:
        for i in set(groups):
            e = a[a.columns[subgroups[sorting]==i]].mean(1)
            e = pd.DataFrame([e for i in range((subgroups[sorting]==i).sum())]).T
            a[a.columns[subgroups[sorting]==i]] = e
    
    fig = sns.clustermap(a, vmin=vmin, vmax=vmax, figsize=size, z_score=0, colors_ratio=0.01, col_cluster=False,col_colors=data, xticklabels=False)
    fig.ax_col_dendrogram.set_visible(False)
    fig.fig.suptitle(title)
    fig.savefig(folder + str(len(set(groups))) + '_clustermap_cobinding_enrichment_andrewplot.pdf')
    plt.show()