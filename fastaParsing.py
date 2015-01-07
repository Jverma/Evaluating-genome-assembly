# -*- coding: utf-8 -*-
#	Parsing a fasta file.
#	Author - Janu Verma
#	jv367@cornell.edu


class FastaParser:
	"""
	Parses a FASTA file to extract the sequences and header information, if any. 

	Parameters
	----------
	fasta_file : Fasta file to be parsed. 


	Example
	-------
	>>> import sys
	>>> input_file = sys.argv[1] 
	>>> out = FastaParser(input_file)
	>>> seqDict = out.sequenceDict()
	>>> print len(seqDict.keys())
	"""

	def __init__(self, fasta_file):
		self.ff = fasta_file
		
	def readFasta(self, fastaFile):
		"""
		Reads and parser the FASTA file. 

		Parameters
		----------
		fastaFile - A FASTA file.

		Returns
		------
		Generator object containing sequences. 
		"""	
		name, seq = None, []
		for line in fastaFile:
			line = line.rstrip()
			if (line.startswith(">")):
				if name: yield (name, ''.join(seq))
				name, seq = line, []
			else:
				seq.append(line)
		if name: yield (name, ''.join(seq))


	def sequenceDict(self):
		"""
		Creates a dictionary of sequences with their header.

		Returns
		-------
		A dictionary of sequences.

		"""
		with open(self.ff) as fastaFile:
			sequences = {}
			for name, seq in self.readFasta(fastaFile):
				sequences[name] = seq
		return sequences		


