## quick hack to reformat vcf files to legend files (read from stdin)
import sys
import gzip
refPanel1kgDir="/research/btc_bioinformatic/operations/Imputation/Data/RefPanel1KG/"

header ='id position a0 a1'# afr.aaf amr.aaf asn.aaf eur.aaf afr.maf amr.maf asn.maf eur.maf' only first four columns essential (at least in impute4)

print(header)
## keeping the info for all snps from 1kg, should use dbSNP/Gnomad, more comprehensive af list
kgpInfo = {line.split()[1]:line for line in open ('%s/ALL_1000G_phase1integrated_v3_chr22_impute_macGT1.legend' % refPanel1kgDir)}
for line in sys.stdin.readlines():
    fields = line.split()
    id, position, ref, alt = [fields[i] for i in [2,1,3,4]]
    if position in kgpInfo:
        print(kgpInfo[position])
    else:
        print(f"{id} {position} {ref} {alt}" + ''.join([' 0.0'] * 8)
    if position in kgpInfo:
        print()

    ## number of missingness should not exceed 50%
    if not line.startswith('#'):
        print Counter(fields[9:162])
        print Counter(fields[162:])
        #print(line, end='')
