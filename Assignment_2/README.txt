  -----------------------------------------------------------------------------------------------------
    This code was created IN COLLABORATION with my colleague Lucía Muñoz Gil.
    Miguel La Iglesia Mirones, Assignment 2, Bioinformatic Programming Challenges.
    Master in Computational Biology, UPM.

    How to run Assignment 2:

    > ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt 
  ---------------------------------------------------------------------------------------------------

This script will make a report named Final_report.txt which contains

    - All networks created after searching protein-protein interactions with the list of coexpressed genes 
    - GO biological processes and KEGG pathways annotations for those networks

Regarding co-expressed gene list from Table S2 of the paper (DOI: 10.1371/journal.pone.0108567) related to Assigment 2, this script was created to have more information
about wether these genes are linked in same regulatory networks or not. We determined if the co-expressed genes are known to bind to one another considering:

    -IntAct data base for searching interactions
    -Uniprot and Kegg for searching protein identifiers and KEGG, GO annotations
    -Filters: quality score (<0.5), physical interaction and detection method (MI:0018)

Results show that out of 168 genes, only 12 genes share the same protein - protein interaction network. Besides that, 150 genes are not included in any network,
that means that no interactions were found for those proteins. Therefore, we suspect that these genes are not actually coexpressed, or at least, they are unlinkely to be in
the same interaction network. Nonentheless, it is true that our results are limited to the data available, after being filtered, in IntAct data base.
