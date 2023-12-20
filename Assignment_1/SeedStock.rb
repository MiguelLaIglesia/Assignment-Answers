=begin
  -----------------------------------------------------------------------------------------------------
    This code was created IN COLLABORATION with my colleague Lucía Muñoz Gil.
    Miguel La Iglesia Mirones, Assignment 1, Bioinformatic Programming Challenges.
    Master in Computational Biology, UPM.
  ---------------------------------------------------------------------------------------------------
=end

# Class SeedStock represents the different seed stocks with their relevant information. 
# It allows to update the grams remmaining and last day planted.

class SeedStock

    # Create an "attribute accessor" (read and write) for those instance properties
    attr_accessor :seed_stock_ID, :mutant_gene_ID, :last_planted, :grams_remaining, :storage
 
    # Predefined values for those instance variables
    def initialize(
        seed_stock_ID: "XXXX", # The ID of the seed stock
        mutant_gene_ID: "ATXXXXX", # The ID of the gene
        last_planted: "dd/mm/yyyy", # Date of the last planted day
        storage:"camaX", # Where it is store
        grams_remaining: 0 # Number of grams in stock
    )

        @seed_stock_ID = seed_stock_ID
        @mutant_gene_ID = mutant_gene_ID
        @last_planted = last_planted
        @storage = storage
        @grams_remaining  = grams_remaining
    end

        # Input: number of grams to be planted
        # Output: number of grams remmaining and last planted day (as instances' properties) are updated.

        def plant_seeds(gr_planted) # Define number of grams planted with gr_planted
            @last_planted = Date.today.strftime("%e/%m/%Y") # Updates last planted day
            if @grams_remaining > gr_planted # Verifies if there is enough grams in stock to be planted
                @grams_remaining -= gr_planted  # Updates grams remaining in stock
            else
                warn "WARNING: we have run out of Seed Stock #{@seed_stock_ID}. #{@grams_remaining} grams were planted" 
                @grams_remaining = 0 # If the grams to be planted is greater than the stock, it warns about it
            end
        end
    
end
