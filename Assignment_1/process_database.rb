=begin
  -----------------------------------------------------------------------------------------------------
    This code was created in collaboration with my colleague Lucía Muñoz Gil.
    Miguel La Iglesia Mirones, Assignment 1, Bioinformatic Programming Challenges.
    Master in Computational Biology, UPM.
  ---------------------------------------------------------------------------------------------------
=end

require 'date' # Library to get today's date
require './SeedStock' # Class SeedStock
require './HybridCross' # Class HybridCross
require './Gene' # Class Gene
require './SeedStockDatabase' # Class SeedStockDatabase

# Gets file names from command line arguments. Verifies a correct number and order of them.

# Verifies number of arguments/ file names
if ARGV.length != 4
    puts "Incorrect number of files passed. There must be 4 files"
    exit 1
end

# Verifies correct file names order
if ARGV[0] == "gene_information.tsv" && ARGV[1] == "seed_stock_data.tsv" && ARGV[2] == "cross_data.tsv" && ARGV[3] == "new_stock_file.tsv"
    gene_info = ARGV[0]
    seed_stock_file = ARGV[1]
    cross_data = ARGV[2]
    new_file_path = ARGV[3]
else 
    puts "Incorrect order of files passed. Usage: ruby process_database.rb  gene_information.tsv  seed_stock_data.tsv  cross_data.tsv  new_stock_file.tsv"
    exit 1
end

# TASK 1. SIMULATION OF 7 GRAMMES SEED PLANTING . RECORDING GENEBANK NEW STATE IN FILE 'seed_stock_data.tsv'

# Create a new instance for class SeedStockDatabase
sesion = SeedStockDatabase.new  

sesion.load_from_file(seed_stock_file)    # Reads seed_stock_data.tsv file to create class SeedStock instances, adding those new instances
                                          # to a hash that allows to access individual SeedStock objects through the stock ID

sesion.stock_instances.each do |_key, gene_object|   # For each class SeedStock instance, it is called a method that updates
    gene_object.plant_seeds(7)                       # the number of grammes remaining and the last date planted
end

sesion.write_database(new_file_path)  # New state of the genebank is written in file 'new_stock_file.tsv' by using properties
                                      # belonging to class SeedStock instances. 

# TASK 2. DETERMINE GENETICALLY LINKED GENES WITH CHI-SQUARE TEST. 

sesion.load_from_file(gene_info)   # Reads 'gene_information.tsv' file to create class Gene instances, adding them
                                   # to a hash that allows to access individual Gene objects based on the gene ID.

# File 'cross_data.tsv' is read to create class HybridCross instances. 
# Then, chi-square test is performed for each of the instances to determine genetically linked genes.
File.open(cross_data, "r").each.with_index do |line, line_num|
    
    next if line_num==0 # Skip file header

    p1, p2, f2_w, f2p1, f2p2, f2p1p2 = line.split()

    crossing = HybridCross.new(parent1_ID:p1, # Instances of class HybridCross are created defining instances'
        parent2_ID:p2,                        # properties by using 'cross_data.tsv' columns
        f2_wildtype:f2_w.to_f,
        f2_p1:f2p1.to_f, # Conversion to float is important to avoid integer rounding when making calculations
        f2_p2:f2p2.to_f,
        f2_p1p2:f2p1p2.to_f)

    x_squared = crossing.chi_squared() # This method is called for each of the instances to calculate
                                       # chi-square value

    if x_squared > 7.82 # If chi-square value is greater than 7.85, then genes are
                        # considered as genetically linked (freedom degrees = 3, p < 0.05)

        # Gets the gene ID from those parents that are genetically linked. It uses the stock ID of parents
        # to get the corresponding instance in class SeedStock, and with that instance it obtains the gene ID.
        linked_gene1 = sesion.get_seed_stock(crossing.parent1_ID).mutant_gene_ID 
        linked_gene2 = sesion.get_seed_stock(crossing.parent2_ID).mutant_gene_ID

        # Gets the gene name from those gene IDs genetically linked. It uses the gene ID to get the corresponding
        # instance in class Gene, and with that instance it obtains the gene name.
        linked_genename1 = sesion.get_gene_info(linked_gene1)  
        linked_genename2 = sesion.get_gene_info(linked_gene2)
        
        # Adds a new property '@linked_to' for those instances in class Gene that have geneticaly linked genes.
        linked_genename1.linked_to = (linked_genename2.name)   
        linked_genename2.linked_to = (linked_genename1.name)
        
        # Prints those linked genes considering the chi-square value 
        puts "Recording: #{linked_genename1.name} is genetically linked to #{linked_genename2.name} with chisquare score #{x_squared}"
    end
end

puts
puts
puts "Final report: \n"
puts

# It parses the instances in class Gene using hash stored in SeedStockDatabase. For those instances
# with genetically linked genes (with atribute ':linked_to'), it prints the corresponding linked gene.

sesion.gene_instances.each do |_key, gene_object|
    if gene_object.respond_to?(:linked_to) && gene_object.linked_to
        puts "#{gene_object.name} is linked to #{gene_object.linked_to}"
    end
end