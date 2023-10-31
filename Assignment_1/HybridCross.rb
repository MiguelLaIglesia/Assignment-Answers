# Class HybridCross represents the different crosses with their relevant information.
# It allows to calculate chi-square in a cross.

class HybridCross

    # Create an "attribute accessor" (read and write) for those instance properties
    attr_accessor :parent1_ID, :parent2_ID, :f2_wildtype, :f2_p1, :f2_p2, :f2_p1p2

    # Predefined values for those instance variables
    def initialize(
        parent1_ID: "XXXX", # The ID of the seed stock of parent 1
        parent2_ID: "XXXX", # The ID of the seed stock of parent 2
        f2_wildtype: 0, # Number of observations for wt phenotype in F2
        f2_p1: 0, # Number of observations for parent 1 phenotype in F2
        f2_p2: 0, # Number of observations for parent 2 phenotype in F2
        f2_p1p2: 0) # Number of observations for parent 1 with parent 2 phenotype in F2

        @parent1_ID = parent1_ID
        @parent2_ID = parent2_ID
        @f2_wildtype = f2_wildtype
        @f2_p1 = f2_p1
        @f2_p2 = f2_p2
        @f2_p1p2  = f2_p1p2
    end

    # No input, it works with intance properties
    # Output: chi-square value ('x_squared')
    # Calculates chi-square value using instance properties that referred to the number of observed F2 phenotypes.
    def chi_squared()
        total_seeds = @f2_wildtype + @f2_p1 + @f2_p2 + @f2_p1p2 # Total number of observations

        # Calculate expected values
        exp_w = total_seeds * 9/16
        exp_het_mutant = total_seeds * 3/16
        exp_hom_mutant = total_seeds * 1/16

        # CalCulate chi-square value
        x_squared = (((@f2_wildtype - exp_w)**2 / exp_w) + ((@f2_p1 - exp_het_mutant)**2 / exp_het_mutant) + ((@f2_p2 - exp_het_mutant)**2 / exp_het_mutant) + ((@f2_p1p2 - exp_hom_mutant)**2 / exp_hom_mutant))
        return x_squared 
    end
end