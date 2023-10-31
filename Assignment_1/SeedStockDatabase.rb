=begin
  -----------------------------------------------------------------------------------------------------
    This code was created IN COLLABORATION with my colleague Lucía Muñoz Gil.
    Miguel La Iglesia Mirones, Assignment 1, Bioinformatic Programming Challenges.
    Master in Computational Biology, UPM.
  ---------------------------------------------------------------------------------------------------
=end

# Class SeedStockDatabase represents the entire data base. It allows to read files to create other classes
# instances. It allows to write updated genebank. It allows to link properties in instances from other classes.

class SeedStockDatabase

    # Create an "attribute accessor" (read and write) for those instance properties
    attr_accessor :stock_instances, :gene_instances

    #Predefined values for those instance variables
    def initialize
        @stock_instances = {}   # This hash allows to access class SeedStock instances through stock ID -> { 'stock ID' : instance }
        @gene_instances = {}    # This hash allows to access class Gene instances through gene ID -> { 'gene ID' : instance }
        @header = "NA"  # Header of file read in 'load_from_file' method
    end

    # Adds class SeedStock instances to hash '@stock_instances'
    def add_seed_stock(seed)
        @stock_instances[seed.seed_stock_ID] = seed
    end

    # Gets class SeedStock instances from hash '@stock_instances' through stock ID
    def get_seed_stock(seed_stock_ID)
        return @stock_instances[seed_stock_ID]
    end

    # Adds class Gene instances to hash '@gene_instances'
    def add_gene(gene)
        @gene_instances[gene.gene_ID] = gene
    end

    # Gets class Gene instances from hash '@gene_instances' through gene ID
    def get_gene_info(gene_ID)
        return @gene_instances[gene_ID]
    end

    # Input A: 'seed_stock_data.tsv' file -> Output: Creates class SeedStock instances, and adds them to the hash that gives access to them through stock ID
    # Input B: 'gene_information.tsv' file -> Output: Creates class Gene instances, and adds them to the hash that gives access to them through gene ID

    def load_from_file(file)
        File.open(file, "r").each.with_index do |line, line_num|
        # Each line in the file will be a class instance
        if file == "seed_stock_data.tsv"
            if line_num==0 # Header is stored in @header, this will be used for method 'write_database'
                @header = line
                next
            end
            # Each column of the file is an instance property value
            seed_ID, gene_id, last_plant_date, storage_place, grams = line.split()
            # Instance is created defining its properties using attr_accessor
            seed_record = SeedStock.new(
                seed_stock_ID:seed_ID,
                mutant_gene_ID:gene_id,
                last_planted:last_plant_date,
                storage:storage_place,
                grams_remaining:grams.to_i)

            # This method is called to add the new class SeedStock instance to the hash stored in SeedStockDatabase object
            add_seed_stock(seed_record)

        elsif file == "gene_information.tsv"

            next if line_num==0
            # Each column of the file is an instance property value
            geneid, gene_name, phenotype = line.split()
            # Instance is created defining its properties using attr_accessor
            gene_record = Gene.new(
                name:gene_name,
                mutant_phenotype:phenotype)
            gene_record.gene_ID=(geneid)

            # This method is called to add the new class Gene instance to the hash stored in SeedStockDatabase object
            add_gene(gene_record)
        end
    end
    end

    # Input: new database file name.
    # Output: Creates and write the updated values of the genebank.

    def write_database(newfile)    
        File.open(newfile, "w") do |file|
            file.puts @header # Header stored when file reading in 'load_from_file'

            # Parse class SeedStock instances using the hash @stock_instances previously created
            # and uses their properties values to update the database file (after simluate planting)
            @stock_instances.each_value do |seed|
                file.puts "#{seed.seed_stock_ID}\t#{seed.mutant_gene_ID}\t#{seed.last_planted}\t#{seed.storage}\t#{seed.grams_remaining}"
            end
        end
    end

end
