import subprocess
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
import pandas as pd
import networkx as nx
import pdb, os

rawData = '/research/gutsybugs/KUMI/Data/RawData/'
flatten = lambda l: [item for sublist in l for item in sublist]

class PlinkInterface:
    def __init__(self, base="merged", wdir='/research/gutsybugs/KUMI/Data/', parseChromoPos=True):
        def extractids(line):
            famID, sampleID = line.split()[:2]
            return (sampleID, famID)
        def extractSNPid(line):
            fields = line.split()
            return (fields[0], int(fields[-1])), fields[1]
        #os.chdir(wdir)
        self.wdir = wdir
        self.base = base
        #infer Family IDs
        self.familyBook = dict([extractids(line) for line in open('%s/%s.fam' % (wdir,base))])
        if parseChromoPos:
            self.chrPos2SNPid = dict(extractSNPid(line) for line in open("%s/%s.map" %(wdir, "alsafar2"))) # alsafar2 -> base # HACK!!
    def readDistanceMatrix(self, upgma=True):
        """
        Performs UPGMA from a distance matrix of samples with genotype-array data
        parses dendrogram into a networkx graph, with which it is easier to identify neighbors
        # expects a plink generated distance matrix: plink --bfile alsafar2 --distance square --out alsafar2
        
        """
        dmfile = "%s%s.dist" % (self.wdir, self.base)
        dmfileIDs = dmfile + ".id"
        self.ids = list(pd.read_csv(dmfileIDs, delimiter='\t', header=None)[1])
        self.dm = pd.read_csv(dmfile, header=None, delimiter='\t')
        self.dm.columns=self.ids
        self.dm = self.dm.set_index(pd.Index(self.ids))
        #print(self.dm)
        if upgma:
            dm = self.dm
            #dms = squareform(dm)
            n = len(self.ids)
            Y = sch.linkage(squareform(self.dm), method='complete') ## consider NJ!
            self.G = nx.DiGraph()
            for index, (a,b,height,clSize) in enumerate(Y):
                ni = n + index
                self.G.add_node(ni, height=height, clSize=clSize)
                self.G.add_edge(ni, int(a)) # weight = height-)
                self.G.add_edge(ni, int(b))
            for i in range(n):
                self.G.node[i]['height'] = 0
                self.G.node[i]['clSize'] = 1
            ## TODO: consider to precalculate allele frequencies for clades of appropriate sizes
    def findClosest(self, nodes, X, increase=30, seedSize=10, minExtent=200):
        """Generator that yields increasing neighborhood supersets for quey nodes:
        Union of kNN for query nodes, with k starting from seedSize, increm. by increase (at least)
        given some nodes that are in our D
        input: nodes - node (sample) names of the query set
        X - the entire training set
        seedSize - 10 Nearest Neighbors 
        minExtend - to avoid to many iterations with very similar neighborhoods, require a minimal increase 

        for each node - get the closest k samples in the training set
        make a union of closest neighbors
        increase k, if the union increased substantially, yield 

        requires readDistanceMatrix to be run first!

        """
        sortedDistances = {node: self.dm[node][X].sort_values() for node in nodes}
        include = seedSize
        lastNeighborsSize = 0
        while True:
            neighbors = set()
            for node in nodes:
                neighbors = neighbors.union(sortedDistances[node][1:include+1].index.values)    
            if len(neighbors) - lastNeighborsSize > minExtent and len(neighbors): 
                lastNeighborsSize = len(neighbors)
                yield neighbors
            if len(neighbors) == len(X):
                break
            include += increase

    def findNeighbors(self, node, nrNeighbors=100, verbose=False):
        node = self.ids.index(node)
        while self.G.node[node]['clSize'] < nrNeighbors:
            if nx.__version__.startswith('1'):
                node = self.G.predecessors(node)[0]
            else:
                ## predecessors returns an iterator in later versions
                node = self.G.predecessors(node).next()
            if verbose:
                print("cluster Size (node %s): %s" % (node, self.G.node[node]['clSize']))
        return [self.ids[leaf] for leaf in self.findLeaves(node)]
    def findNeighborSets(self, node, nrNeighborsMin=100, nrNeighborsMax=None):
        node = self.ids.index(node)
        while True: #self.G.node[node]['clSize'] < nrNeighborsMax:
            if not self.G.predecessors(node):
                break
            else:
                node = self.G.predecessors(node)[0]
                if self.G.node[node]['clSize'] > nrNeighborsMin:
                    yield [self.ids[leaf] for leaf in self.findLeaves(node)]  
    def findNeighborSetsMulti(self, nodes, nrNeighborsMin=100, nrNeighborsMax=None):
        """TODO: find suitable neighbors for a set of samples, not just one"""

    def findLeaves(self, node):
        successors = list(self.G.successors(node))
        if len(successors) == 0:
            return [node]
        return flatten([self.findLeaves(successor) for successor in successors])

    def extractVariants(self, variantFile, outfile):
        cmd = 'plink --bfile %s%s --extract %s --make-bed --out %s' % (self.wdir, self.base, variantFile, outfile)
        print(cmd)
        subprocess.call(cmd, shell=True)
    def keepNeighbors(self, sample, base=None):
        if base is None: base == self.wdir + self.base
        cmd = 'plink --bfile %s --keep %s --make-bed --out %s' % (base, sample.neighborsFile, sample.variantOverlapNeighbors)
        print(cmd)
        subprocess.call(cmd, shell=True)

    def calculateAlleleFreqs(self, variantFile):
        cmd = 'plink --bfile %s --freq --out %s' % (variantFile, variantFile)
        print(cmd)
        subprocess.call(cmd, shell=True)
        return pd.read_csv('%s.frq' % variantFile, sep=r'\s*', index_col=1)

class Individual:
    def __init__(self, id, vcf):
        self.id = id
        self.vcf = vcf
        self.dir = rawData + id
        self.variantFile = '/tmp/Rare/variants_%s.txt' % id
        self.neighborsFile = '/tmp/Rare/neighbors_%s.txt' % id  #list of sample ids
        self.variantOverlap = '/tmp/Rare/variantsOverlap_%s' % id #bed/bim/fam + frq
        self.variantOverlapNeighbors = '/tmp/Rare/variantsOverlapNeighbors_%s' % id #bed/bim/fam + frq

    def readVariants(self):
        def processline(line):
            if line.startswith('#'): return None, False
            chrom, pos, rsid, ref, alt, qual, filter, info = line.split()[:8]
            return (chrom.lstrip('chr'), int(pos)), filter == 'PASS' and rsid.startswith('rs')

        A0 = [processline(line) for line in open('%s/%s'%(self.dir, self.vcf))]
        A =  [variant for (variant, ok) in A0 if ok]
        self.variants = set(A)

    def writeVariants(self, plink):
        with open(self.variantFile, "w") as vf:
            for variant in self.specific: ## can be extended to all!!
                if variant in plink.chrPos2SNPid:
                    vf.write('%s\n' % plink.chrPos2SNPid[variant]) # kgp8671108 instead of chrom+pos
        print("Wrote %s" % self.variantFile)

    def writeNeighbors(self, plink):
        neighbors = plink.findNeighbors(int(self.id), 50)  # sample.id) ## TODO: FIX!!!
        with open(self.neighborsFile, "w") as sf:
            for neighbor in neighbors:
                sf.write("%s\t%s\n" % (plink.familyBook[str(neighbor)], neighbor))

    def calculateAlleleFreqs(self, rareThreshold = 0.05):
        subpopulationAFs = plink.calculateAlleleFreqs(self.variantOverlapNeighbors)
        populationAFs = plink.calculateAlleleFreqs(self.variantOverlap)
        populationAFs.columns = [col+'2' for col in populationAFs.columns.values]
        mafTable = pd.concat([subpopulationAFs,populationAFs], axis=1, join='inner')
        mafTable['DIFF'] = mafTable.apply(lambda row: row.MAF - row.MAF2, axis=1)
        mafTable[(mafTable['MAF']<rareThreshold) & (mafTable['MAF2']>rareThreshold)]
        mafTable.sort('DIFF', inplace=True)

    def individualVariants(self, others):
        """given other variants"""
        self.specific = self.variants - others

    def selectKin(self, k=100):
        """use AESA and IBS distance matrix to infer nearest k samples)"""
        self.subpop = None

class Cohort:
    def __init__(self, ids, vcffile='_recalibrated_variants.vcf'):
        self.individuals = [Individual(str(id), '%s%s'%(id, vcffile)) for id in ids]

    def readVariants(self):
        for sample in self.individuals:
            sample.readVariants()

    def calcIndividualVariants(self):
        for i in range(len(self.individuals)):
            s = self.individuals.pop(0)
            s.individualVariants(set().union(*[c.variants for c in self.individuals]))
            self.individuals.append(s)

    def writeVariants(self, plink):
        for sample in self.individuals:
            sample.writeVariants(plink)

    def extractVariantOverlap(self, plink):
        for sample in self.individuals:
            plink.extractVariants(sample.variantFile, sample.variantOverlap)

    def keepNeighbors(self, plink):
        for sample in self.individuals:
            sample.writeNeighbors(plink)
            plink.keepNeighbors(sample, base=sample.variantOverlap)
    def calculateAlleleFreqs(self):
        for sample in self.individuals:
            sample.calculateAlleleFreqs()

## Please name vcf files consistently, than this can be more elegant in a loop here:
if __name__=="__main__":
    cohort = Cohort([10745,12584,13076])

    #See  plink filter input keep, extract
    #plink --bfile mydata3 --extract plink.prune.in --make-bed --out mydata4
    #plink --bfile data --keep mylist.txt --make-bed --out data_keep

    #plink = PlinkInterface(base='alsafar2', wdir='/home/ahenschel/Dropbox/KU-MI 20170504 Anthropology study/Data')
    plink = PlinkInterface(parseChromoPos=False) ## going with default while on campus
    #plink.readDistanceMatrix()

    #cohort.readVariants()
    #cohort.calcIndividualVariants()
    #cohort.writeVariants(plink)
    #cohort.extractVariantOverlap(plink)
    #cohort.keepNeighbors(plink)
    cohort.calculateAlleleFreqs()

