# -*- coding: utf-8 -*-
#	Computing statistics of Assembled Contigs.
#	Author - Janu Verma
#	jv367@cornell.edu



import sys
import numpy as np 
import csv
import matplotlib 
import pylab
from fastaParsing import FastaParser


class AssemblyStatistics:
	"""
	Computes basic statistics of the assembled contigs. 

	Parameter
	---------
	inputFile : The FASTA file containing contigs. 

	Example
	-------
	>>> import sys
	>>> import numpy as np
	>>> inputFile = sys.argv[1]
	>>> out = AssemblyStatistics(inputFile) 
	>>> N50 = out.N50()
	>>> minContigLength = out.minContigLength()
	>>> maxContigLength = out.maxContigLength()

	"""
	def __init__(self, inputFile):
		self.inFile = inputFile
		self.fastaInfo = FastaParser(self.inFile)
		self.contigsInfo = self.fastaInfo.sequenceDict()


	def scores(self):
		"""
		Compute the basic statistics. 

		Returns
		-------
		Dictionary of the basic statistics of the assembly.
		"""
		seqLengths = []
		for x in self.contigsInfo.keys():
			seq = self.contigsInfo[x]
			seqLengths.append(len(seq))

		seqLengths = sorted(seqLengths)	
		max_length = max(seqLengths)
		min_length = min(seqLengths)
		mean_length = np.mean(seqLengths)	


		midLength = sum(seqLengths)/2

		computedMidLength = 0
		l50 = 0
		n50 = 0
		for i,x in enumerate(seqLengths):
			if (midLength < computedMidLength):
				n50 = i
				l50 = x 
				break
			computedMidLength += x

		scoresDict = {'number_of_contigs':len(seqLengths), 'smallestContig':min_length, 'meanContig':mean_length, 
		'n50':n50, 'l50':l50, 'largestContig':max_length, 'lengthOfAssembly':sum(seqLengths)}
		return scoresDict


	def N50(self):
		"""
		Computes the N50 score of the assembly. 

		Returns
		-------
		Numerical value of the N50 statistics.
		"""
		stats = self.scores()
		return stats['n50']

	def L50(self):
		"""
		Computes the L50 score of the assembly. 

		Returns
		-------
		Numerical value of the L50 statistics.
		"""
		stats = self.scores()
		return stats['l50']


	def maxContigLength(self):
		"""
		Computes the length of the largest contig in the assembly. 

		Returns
		-------
		Numerical value of the max length.
		"""
		stats = self.scores()
		return stats['largestContig']

	
	def minContigLength(self):
		"""
		Computes the length of the smallest contig in the assembly. 

		Returns
		-------
		Numerical value of the min length.
		"""
		stats = self.scores()
		return stats['smallestContig']


	def meanContigLength(self):
		"""
		Computes the mean length of the contigs in the assembly. 

		Returns
		-------
		Numerical value of the mean length.
		"""
		stats = self.scores()
		return stats['meanContig']



	def lengths(self):
		"""
		Write the contig lengths in a csv file. 

		Returns
		-------
		A csv file (contigLengths.csv) containing the lengths of the contigs.
		"""
		output_file = csv.writer(open('contigLengths.csv', 'wb'))
		for i,x in enumerate(self.contigsInfo.keys()):
			seq = self.contigsInfo[x]
			l = len(seq)
			output_file.writerow([i,l])



	def histogramOfContigLengths(self):
		"""
		Plots a histogram of the contig lengths. 

		Returns
		------
		A histogram of contig sizes, saved in the file contig_histogram.png
		"""
		seqLengths = []
		for x in self.contigsInfo.keys():
			seq = self.contigsInfo[x]
			seqLengths.append(len(seq))

		seqLengths = sorted(seqLengths)
		l = seqLengths[0:540000]
		matplotlib.use('Agg')
		pylab.hist(l, bins=50)
		pylab.title("Contigs historgram")
		pylab.xlabel('Sequence Length (bp)')
		pylab.ylabel('Count')
		pylab.savefig('contig_histogram.png')	


	def boxPlot(self):
		"""
		Plots a box-plot of the contig lengths. 

		Returns
		------
		Box plot of contig sizes, saved in the file contig_boxplot.png
		"""
		seqLengths = []
		for x in self.contigsInfo.keys():
			seq = self.contigsInfo[x]
			seqLengths.append(len(seq))

		pylab.boxplot(seqLengths)
		pylab.savefig('contig_boxplot.png')	


	def assemblyCoverage(self, genomeSize):
		"""
		Computes the coverage of the assembly. 

		Parameter
		---------
		genomeSize : (Expected )Size (in bp) of the genome. 

		Returns
		-------
		Numerical value of the fraction of the genome covered by the assembly. 
		"""
		stats = self.scores()
		assemblySize = stats['lengthOfAssembly']
		return float(assemblySize)/genomeSize


