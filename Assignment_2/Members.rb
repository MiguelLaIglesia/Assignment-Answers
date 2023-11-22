# == Members
#
# This is a representation of proteins
# with a Uniprot ID, gene ID and its direct interactors
#
# == Summary
# 
# This can be used to represent members that participate in a protein-protein interaction network
#

class Members

    # Get/Set the protein's ID
    # @!attribute [rw]
    # @return [String] The Uniprot ID

    attr_accessor :uniprot_id

    # Get/Set the proteins' direct interactors
    # @!attribute [rw]
    # @return [array<Members>] The array of direct interactors
    
    attr_accessor :direct_interactors

    @@coexpresed_members= []
    @@all_members = Hash.new
    @@total_num_members = 0

    # Create a new instance of Members
  
    # @param uniprot_id [String] the Uniprot ID of the member as a String
    # @return [Members] an instance of Members

    def initialize( uniprot_id: "IDXXX" )
        @uniprot_id = uniprot_id
        @direct_interactors = []                # direct_interactors [Array<Members>] the direct interactors of the member as an array of Members' instances
        @@all_members[@uniprot_id] = self       # all_members [Hash] the totality of the members accessed by the uniprot_ID
        @@total_num_members +=1                 # total_num_members [Integer] the total number of members
    end

    # Set the gene_id of the member
    # @param gene_name [String] the gene name of the member
    # @return [void]

    def gene_id=(gene_name)
        @gene_id = gene_name
    end

    # @!method no_params_method
    # Get the gene_id of the member
    # @return [String] the gene ID of the member

    def gene_id
        @gene_id
    end

    # @!method no_params_method
    # Get all members created with coexpressed gene list
    # @return [Array<Members>] the list of coexpressed members

    def self.all_coexpresed_members
        @@coexpresed_members
    end

    # @!method no_params_method
    # Get all members created
    # @return [Array<Members>] the list of coexpressed members
    def self.all_members
        @@all_members
    end

    # @!method no_params_method
    # Get total number of members
    # @return [Integer] the number of members
    def self.number_of_members
        @@total_num_members
    end

    # Creates members while reading each line of a file, and save them in a class atribute
    # @param filename [String] the name of the file to read gene identifiers
    # @return [void]

    def self.read_from_file(filename)
        coexpressed_file = File.open(filename, 'r')
        coexpressed_file.readlines.each do |line|
            locus_name=line.chomp

            if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Regex for Arabidopsis thaliana genes
                abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
            end

            result = togo_search("uniprot", locus_name,"/accessions") # Search that gene in togo REST API
            if result.is_a?(Array) && result.any?
                uniprot_id = result.first.first  # First uniprot ID is stored 
            else
                puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
                next
            end
            member = AnnotatedMembers.new(uniprot_id: uniprot_id) # Create AnnotatedMember for that gene
            member.gene_id=(locus_name)
            @@coexpresed_members << member
        end
    end

    
    
    # Finds and stores protein interactors for a given protein by accessing IntAct database
    # @param intact_address [String] the IntAct URL
    # @param species [String] the specie of the organism to search
    # @param formato [String] the output file format retrieved in the search 
    # @return [void]

    def find_interactors(intact_address=INTACT_BASE_ADDRESS, species=SPECIES, formato=TAB25_FORMAT)
        intact_address = "#{intact_address}search/interactor/#{@uniprot_id}/?query=#{species}&format=#{formato}" # Access through PSICQUIC REST API
        response = rest_api_request(intact_address)
        if response.empty? 
          @direct_interactors = "Response Not Available in IntAct"
          return 1
        end
        response.body.each_line do |line|
            values = line.chomp.split("\t")

            # Filtering by confidence score
            instact_miscore = extract_xref(values[14]).to_f
            next if instact_miscore < 0.5
            
            # Filtering by quality of/trust in tecnology
            interaction_detection_method = extract_xref(values[6])
            next if interaction_detection_method == "MI:0018"

            # Filtering by type of interaction
            type_int = extract_xref(values[11])
            next if type_int != "MI:0915"
            
            # Select the interactor that is not the query and store it in direct interactors
            [0,1].each do |id| 
                interactor = extract_xref(values[id])
                if !interactor.nil? && !interactor.include?(@uniprot_id) && interactor.match(/[OPQ][0-9][A-Z0-9]{3}[0-9]$/) # Regex for uniprot ID
                    if @@all_members.key?(interactor) # Asking if that member already exists
                        @direct_interactors << @@all_members[interactor]
                    else
                        interactor = AnnotatedMembers.new(uniprot_id: interactor) 
                        @direct_interactors << interactor
                    end 
                end
            end
        end
    end

end