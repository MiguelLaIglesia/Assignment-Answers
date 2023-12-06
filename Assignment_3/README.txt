  -----------------------------------------------------------------------------------------------------
    This code was created IN COLLABORATION with my colleague Lucía Muñoz Gil.
    Miguel La Iglesia Mirones, Assignment 3, Bioinformatic Programming Challenges.
    Master in Computational Biology, UPM.

    How to run Assignment 3:

    > ruby main.rb ArabidopsisSubNetwork_GeneList.txt
  ---------------------------------------------------------------------------------------------------

This script will take Arabidopsis thaliana gene names from list and will search embl files in order to find all repetetive 'CTTCTT' regions inside exons.
In this work, several 'CTTCTT' regions together (one after the other) were considered as different (and not only one) repetitive regions. We took this decision
regarding biologists want inserts to go in each of those regios, anyway, it would be interesting a further discussion with biologists to know how they require it.