How to run Assignment 2 for Bioinformatic Programming Challenges:

> ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt 

This script will make a report named Final_report.txt which contains
    - All networks created after searching protein-protein interactions with the list of coexpressed genes 
    - GO biological processes and KEGG pathway annotations for those networks

Regarding the co-expressed gene list from Table S2 of the paper (DOI: 10.1371/journal.pone.0108567) related to Assigment 2, this script was run to have some clues
about the linking of genes in regulatory networks. Therefore, we determined if the co-expressed genes are known to bind to one another considering:
    -IntAct data Base for searching interactions
    -Uniprot and kegg for searching protein identifiers and KEGG,GO annotations
    -Filters: quality score (<0.5), physical interaction and detection method (MI:0018)

Results show that out of 168 genes, are only 12 genes share the same protein - protein interaction network. Besides that, 150 genes are not included in any network,
that means that no interactions were found for those proteins. Therefore, we suspect these genes are actually coexpressed, or at least, they are unlinkeli to take in
the same interaction network. Nonentheless, it is true that our results are limited to the data available, after being filtered, in IntAct data base.

