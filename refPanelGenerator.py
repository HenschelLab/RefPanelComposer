"""
Plink interface to get distance matrix of all samples (both target and ref)
For starters, only do chromosome 22, Proof of Concept

Panels are generated in impute format by extracting columns from the "merged" reference panel,
that includes UAE (u), Qatar (q) and 1kGenome (k) -> uqk_chr22.hap

Algorithm: 

Data Preparation:
for UAE, Q, do phasing (shapeit), select biallelic variants only
merge vcf files of the ref panels to obtain UQK
convert to bed/bim/fam (plink)
Distance Matrix calculation for all (!) samples. I use plink/1-IBS on 1-IBS

given one or more query genotype-arrays
split UAE data set into training and test

create reference panel from a subset of UAE-training + Qatar + 1KG, using findClosest.
zcat uqk_chr22.hap.gz | cut -d' ' -f1-16,18-52,54-65,67-79,82-127 ... | gzip > uqk_chr22_00_3283.hap.gz
(where 00 is the cross validation fold and 3282 is the size of the reference panel, here the total)

## Run like this on HPC:
source activate bionew
Out[24]: cd '/research/btc_bioinformatic/operations/Imputation/Data/Combined'
In [25]: run ../../refPanelGenerator.py 
(or fix path for some of the Imputation files)
----------------------------------------------------
class PlinkInterface:
    def findClosest(genotype-arrays):     
    see more in PlinkInterface.findClosest docstring

"""
import pdb
import gzip
from PlinkInterface import PlinkInterface, Cohort, Individual
import allel
from sklearn.model_selection import ShuffleSplit
import subprocess
from filterMAC import filterMAC
computeFacility = ['HPC', "BTC"][0]

if computeFacility =='HPC':
    basedir = "/research/btc_bioinformatic/operations/Imputation/Data"
    datadir=f"{basedir}/Combined/"
    refPanel1kgDir=f"{basedir}/RefPanel1KG/"
    geneticMapDir=f"{basedir}/GeneticMaps/"
    targetDir = f"{basedir}/UAE"
    #genericRef = "uqk_arabSNPs_chr22.hap.gz"

def rangifyer(i, offset=0):
    import itertools
    for a, b in itertools.groupby(enumerate(i), lambda x: x[1] - x[0]):
        b = list(b)
        if b[0][1] == b[-1][1]:
            yield str(b[0][1] + offset)
        else: 
            yield "%s-%s" % (b[0][1] + offset, b[-1][1] + offset)

class Imputation:
    def __init__(self, X, panelBase, fold, phasedTargets='phased_chr%s.shapeit', chromosome=22):
        """Expecting as input phased target (query) samples (option known_haps_g) """
        phasedTargets=phasedTargets % chromosome
        self.chromosome = chromosome
        self.genericRef = panelBase+'.hap.gz'
        self.vcfAll = panelBase+'.vcf.gz'
        ## Setting names
        self.fold = fold
        panel0 = '%s_%02d_%03d' % (panelBase, fold, len(X)) ##  note fold and size of training set into filename
        self.panel = f'{panel0}.hap.gz'
        self.legend = f'{panel0}.legend.gz'
        self.panel1 = f'{panel0}_MAC%s.hap.gz'
        self.legend1 = f'{panel0}_MAC%s.legend.gz'
        self.phasedTargetSamples = f'{targetDir}/{phasedTargets}.haps'
        sampleIds = f'{targetDir}/{phasedTargets}.sample'
        ## only needed for checking that samples in haps file are in same order as combined data 
        self.phasedTargetSampleIds = [line.strip().split()[0] for line in open(sampleIds)][2:]        
        self.knownHaps = '%sFold%02d.haps' % (self.phasedTargetSamples[:-4], self.fold)
        self.knownHapsMasked = '%s_MSKD.haps' % (self.knownHaps[:-5])
        


    def filterRefLowMAC(self, minmac=5):
        ## remove snps with low Minor Allele Count (<minmac)
        ## Careful: requires the deletion of the corresponding line in the legend file, I guess
        ## simultaneous sweep through both files (panel/legend)
        filterMAC(self.panel, self.legend, self.panel1%minmac, self.legend1%minmac, minmac)

    def vcf2legend(self):
        """"Produces complete legend for all variants in initial vcf."""
        import os
        if os.path.exists(self.legend):
            print("Warning, legend file exists! Skipping.")
            return
        with gzip.open(self.legend, 'wt') as out:
            snpCounter = 0
            for line in gzip.open(self.vcfAll):
                line = line.decode()
                if line.startswith('#'): continue
                fields = line.strip(). split()
                rsid, position, ref, alt = [fields[i] for i in [2,1,3,4]]
                if not rsid.startswith('rs'):
                    rsid = 'arab%06d' % snpCounter
                    snpCounter += 1
                print(f"{rsid} {position} {ref} {alt}", file=out)

    def doubleCheckSampleOrder(self, targetIDs):
        if targetIDs != self.phasedTargetSampleIds:
            pdb.set_trace()
    def writeTestset(self, yIdx, maskBimFile='ABchr'):
        """creates the test set (y, aka target samples) that is to be imputed as hap file. They are cut out from the file that contains all (phased) UAE
        haplotypes
        """
        if maskBimFile: ## taking mask from genotype array bim file -> all genotyped positions
            mask = [int(line.split()[3]) for line in open(f'{targetDir}/Bim/{maskBimFile}{self.chromosome}.bim')]
            with open(self.knownHapsMasked, 'wt') as khm:
                for line in open(self.phasedTargetSamples):
                    pdb.set_trace()

        ## cut from UAE samples, always include first 5 columns 
        cmd = "cut -d' ' -f 1,2,3,4,5,%s %s > %s" % (','.join(rangifyer(yIdx, offset=5)), \
                                    self.phasedTargetSamples, self.knownHaps)


    def generateRefPanel(self, X, allSamples, dryRun=False):
        """produce .hap files that are used as ref panels in impute2
        Alternatively look at impute's merge option, takes much of the file format headache and 
        does a mutual impute on the refpanels        
        """        
        xIdx = sorted([allSamples.index(x)+1 for x in X])
        ## select columns from all-including hap file using linux' cut
        cmd = "zcat %s | cut -d' ' -f" % (self.genericRef) + ','.join(rangifyer(xIdx)) + '| gzip > %s' % (self.panel)
        print(cmd)
        ## select columns of phased target samples (y)
        if not dryRun:
            subprocess.call(cmd, shell=True)
                
    def runImpute(self, chrom=22, regionStart=16000000,  regionEnd =21000000):
        mapFile = "genetic_map_3col_GRCh37_chr%s.txt" % chrom
        geneticMap = "%s/%s" % (geneticMapDir, mapFile)
        prephased = "%s" % (combinedDir, prephasedTrain)
        outfile = "%s_%s_CV%s" % () ## TODO!!! 
        cmd = f"impute2 -m {geneticMap} -known_haps_g {prephased} -h {self.panel} -l {self.legend} -Ne 20000 -int {regionStart} {regionEnd} -o {outfile} -allow_large_regions"
        print(cmd)
        #subprocess.call(cmd)
        
chrom = 22

uaeCount = 153
qtrCount = 1005
kgpCount = 2504

panelBase = 'uqk_all_chr22' # 'uqk_arabSNPs_chr22'
vcf=f"{panelBase}.vcf.gz" ## to be updated with merged file

## get SampleIds FullGenome(sampleFile):
data = allel.read_vcf_headers(f"{datadir}/{vcf}")

## TODO: filter dataset by SNPs from one chromosome, for starters
pi = PlinkInterface(wdir=datadir, base=panelBase, parseChromoPos=False)
pi.readDistanceMatrix(upgma=False)

## increasing neighbor set
splitter = ShuffleSplit(n_splits=5,  test_size=.2, random_state=0)
#splitter.get_n_splits(data.samples)
fold = 0
## only split UAE samples
## Cross validation, each round leaves some UAE test samples out
for trainIdx0, testIdx in splitter.split(data.samples[:uaeCount]):
    ## choose around 30 UAE samples as test set, remaining 120 + others as superset from which training sets are selected using 
    trainIdx = list(trainIdx0) + list(range(len(data.samples[uaeCount:])))
    Xsuper = [sample for idx, sample in enumerate(data.samples) if idx in trainIdx]
    y = [sample for idx, sample in enumerate(data.samples) if idx in testIdx]

    firstTime=True
    for X in pi.findClosest(y, Xsuper):
        impute = Imputation(X, panelBase, fold, chromosome=22)
        #impute.doubleCheckSampleOrder(data.samples[:uaeCount])
        if firstTime:
            impute.vcf2legend() ## complete legend of all 
            impute.writeTestset(testIdx) ## shouldn't be necessary for every neighborhood, only once, well...
            firstTime = False
                
        impute.generateRefPanel(X, data.samples, dryRun=True)
        impute.filterRefLowMAC()
        break # just try with a small X
    break
    fold += 1

"""
[ahenschel@sub2 UAE]$ zcat phased_chr22.shapeit.vcf.gz |grep -v '^#' |wc -l
242435 (nr of variants in chr 22)
[ahenschel@sub2 UAE]$ zcat gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_onlySNP_chr22.vcf.gz |grep -v '^#' |wc -l
246403
[ahenschel@sub2 UAE]$ #zcat gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_onlySNP_chr22.vcf.gz |grep -v '^#' |wc -l
[ahenschel@sub2 UAE]$ zcat gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_onlySNP_chr22.vcf.gz |grep -v '^#' |cut -d$'\t' -f3 |grep -v 'rs'
[ahenschel@sub2 UAE]$ zcat gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_onlySNP_chr22.vcf.gz |grep -v '^#' |cut -d$'\t' -f3 |head -3
rs62224609
rs3013006
rs372030108


Notes:
## file prep (chronologically)
[ahenschel@sub2 UAE]$ ls -lt *.vcf.gz

---------------
## splitting (with tabix?, filter: SNPs only)
-rw-r--r-- 1 ahenschel henschel-lab  204459997 Aug 26 18:48 gathered_vcfs_recal_snp_indel_snpeff_dbsnp_clinvar_onlySNP_chr22.vcf.gz
--------------
## From the SNPs, just choose biallelic
-rw-r--r-- 1 ahenschel henschel-lab 1090119532 Aug 27 12:54 biallelic_chr1.vcf.gz ...
-rw-r--r-- 1 ahenschel henschel-lab  194888614 Aug 27 12:49 biallelic_chr22.vcf.gz
---------------
## Phasing with shapeit, now files are in operation/Imputation/Data/UAE 
## see operation/Phasing/shapeitCmd.txt
-rw-r--r-- 1 ahenschel henschel-lab   29771813 Aug 30 15:47 phased_chr10.shapeit.vcf.gz ...
-rw-r--r-- 1 ahenschel henschel-lab   25550486 Aug 30 15:50 phased_chr9.shapeit.vcf.gz

# vcf to ped
plink --vcf uqk_chr22.vcf.gz --recode --out uqk_chr22

# ped -> bed
# plink --file uqk_arabSNPs_chr22 --make-bed --out uqk_arabSNPs_chr22
 
## this refpanel is to be improved, but have to run the original as well to have a baseline
refpanel = "%s/ALL_1000G_phase1integrated_v3_chr%s_impute_macGT1.hap.gz" % (refPanel1kgDir, chrom)
refpanelLegend = "%s/ALL_1000G_phase1integrated_v3_chr%s_impute_macGT1.legend.gz" % (refPanel1kgDir, chrom)
mapFile = "genetic_map_3col_GRCh37_chr%s.txt" % chrom
geneticMap = "%s/%s" % (geneticMapDir, mapFile)
#prephased = "%s" % (combinedDir, prephasedTrain)
#outfile = "%s_%s_CV%s" % () ## TODO!!!

#f"impute2 -m {geneticMap} -known_haps_g {prephased} -h {refpanel} -l {refpanelLegend} -Ne 20000 -int {regionStart} {regionEnd} -o {outfile} -allow_large_regions"

## Ref Panel format (to be emulated from combined vcf)
ahenschel@sub3 RefPanel1KG]$ zcat ALL_1000G_phase1integrated_v3_chr22_impute_macGT1.hap.gz |head -1
0 0 0 1 0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ... (2184)
394274 lines (SNPs) not counting 1 header line

zmore ALL_1000G_phase1integrated_v3_chr22_impute_macGT1.legend.gz |head -3
id position a0 a1 afr.aaf amr.aaf asn.aaf eur.aaf afr.maf amr.maf asn.maf eur.maf
rs149201999 16050408 T C 0.0955 0.0497 0.0437 0.058 0.0955 0.0497 0.0437 0.058
394275 lines (SNPs) with one 1 header line


nohup time impute2 -m /bmshare/ahenschel/References/Shapeit/genetic_map_3col_GRCh37_chr22.txt \
        -known_haps_g ABchr22.shapeit.haps\
        -h ALL_1000G_phase1integrated_v3_chr22_impute_macGT1.hap.gz\
        -l ALL_1000G_phase1integrated_v3_chr22_impute_macGT1.legend.gz\
        -Ne 20000\
        -int 16000000 19000000\
        -o ABchr22.pos16M-19M.impute2\
        -allow_large_regions\
        -seed 367946 > nohup_impute2_chr22.out &
        
        
        
        
        
Aux:
            ## Creating Legend file
            with open("%s/%s" %(datadir, self.legend), "w") as legfile:
                for xidx in xIdx:
                    print(allSamples[xidx], file=legfile)

        
        
        
        
        """