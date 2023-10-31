# Class Gene represents the different genes with their relevant information.
# It allows to verify proper gene ID format and to store genetically linked genes.

class Gene

    # Create an "attribute accessor" (read and write) for those instance properties
    attr_accessor :name, :mutant_phenotype
    
    #Predefined values for those instance variables
    def initialize(
        name: "XXXX", # Gene name
        mutant_phenotype: "NA" # Description of the phenotype
    )
        @name = name
        @mutant_phenotype = mutant_phenotype
    end

    # Gets 'gene_ID' property
    def gene_ID
        @gene_ID
    end

    # Set 'gene_ID' property to the value passed to the method, verifying a proper format of the value
    def gene_ID=(string)
        unless string.match(/A[Tt]\d[Gg]\d\d\d\d\d/) # .match allows to validate that regular expression
            warn "CODE STOP: Gene Identifier format is incorrect. It should have the format '/A[Tt]\d[Gg]\d\d\d\d\d/' "
            exit 1
        end
        @gene_ID = string # Set gene_ID instance property to the value passed to the method
    end

    # Gets 'linked_to' property
    def linked_to
        @linked_to
    end

    # Set 'linked_to' property to the value passed to the method
    def linked_to=(new_gen)
        if @linked_to # This is only in case the gene was genetically linked to several genes
            @linked_to += ",#{new_gen}"
        else
            @linked_to = new_gen # Sets a value to the property
        end
    end

end