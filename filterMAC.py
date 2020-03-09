import gzip
def filterMAC(origPanel, origLegend, newPanelFile, newLegendFile, minmac=5):
    with gzip.open(newLegendFile, 'wt') as newLegend:
        print('id position a0 a1', file=newLegend)
        with gzip.open(newPanelFile, 'wt') as newPanel:
            with gzip.open(origPanel, 'rt') as panel:
                for lline in gzip.open(origLegend, 'rt'):
                    ## this will crash if panel has fewer lines than legend! Shouldn't
                    panelLine = panel.readline()
                    mac = sum(map(int, panelLine.rstrip().split()))
                    if mac >= minmac:
                        print(lline, end='', file=newLegend)
                        print(panelLine, end='', file=newPanel)

if __name__=="__main__":
    ## run in /research/btc_bioinformatic/operations/Imputation/Data/Combined with
    ## run ../../filterMAC.py
    filterMAC('uqk_all_chr22_00_235.hap.gz', 'uqk_all_chr22_00_235.legend.gz',\
              'uqk_all_chr22_00_235_TEST.hap.gz', 'uqk_all_chr22_00_235_TEST.legend.gz')
    
