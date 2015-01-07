# Evaluating-genome-assembly
computation of some statistics used to evaluate a genome assembly. 

<b>Dependencies:</b>
        numpy
        matplotlib


<b>Usage:</b>

        import sys
        from matplotlib import pyplot
        from stats import AssemblyStatistics
        
        # the input contig file, in FASTA format. 
        inputFile = sys.argv[1]
        
        
        out = AssemblyStatistics(inputFile)
        
        # L50 of the assembly
        l50 = out.L50()
        
        # N50 of the assembly
        n50 = out.N50()
        
        # size of the largest contig
        largestContig = out.maxContigLength()
        
        # size of the samallest contig
        smallestContig = out.minContigLength()
        
        # mean contig size 
        meanContig = out.meanContigLength()
        
        # genome coverage
        coverage = out.assemblyCoverage()
        
        # histogram of the contig lengths
        out.histogramOfContigLengths()
        
        # box plot of contig lengths
        out.boxplot()
        
