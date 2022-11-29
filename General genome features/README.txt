##########################################################################################################
#############################                                               #############################
#############################    README file for General genome features    #############################
#############################                                               #############################
##########################################################################################################

**BUSCO - Benchmarking sets of Universal Single-Copy Orthologs**
###___command for BUSCO__### run_BUSCO.py -i input.fas -o output.txt file -l fungi_odb9 -m proteins -c 40
##For BUSCO runs fasta files for each species should be provided as a seperate file
##used proteome fasta in this case, for which outputs for each species will be saved in their own directory
##the output file we need, will have ".txt" extension (e.g. "short_summary_inputfungi.txt")
##to extract the "short_summary___.txt" files from multiple directories into a single directory, use "busco_outmerg.sh"

**Concatenated supermatrix for inferring species tree**
##Scripts used for concatenation of trimmed alignments into supermatrix and creating the partition file for species tree inference was taken from Zsolt Merenyi (https://github.com/zsmerenyi/)

**COMPARE - COMparative Phylogenomic Analysis of tRait Evolution **
(Nagy et al., 2014 - https://www.nature.com/articles/ncomms5471)
##scripts for COMPARE pipeline were based on https://github.com/laszlognagy/COMPARE and https://github.com/zsmerenyi/compaRe

