# == AnnotatedMembers
#
# This is a representation of members functionally annotated
# with a GO or KEGG annotations, but many others could be included by just defining the proper method
#
# == Summary
# 
# This can be used to represent members that participate networks considering their functional annotations
#


class AnnotatedMembers < Members

    # Get/Set the member GO annotations
    # @!attribute [rw]
    # @return [Hash] The GO ID and GO term
    attr_accessor :go_IDs_terms

    # Get/Set the KEGG's gene ID
    # @!attribute [rw]
    # @return [String] The gene ID in KEGG
    attr_accessor :kegg_gene

    # Get/Set the member KEGG annotations
    # @!attribute [rw]
    # @return [Hash] The KEGG ID and KEGG term
    attr_accessor :kegg_ID_pathway


    # Create a new instance of AnnotatedMembers    
    # @return [AnnotatedMembers] an instance of AnnotatedMembers
    def initialize(**args)
        super
        @go_IDs_terms = {}         # go_IDs_terms [Hash] the GO IDs and their corresponding description/term
        @kegg_gene = ""            # kegg_gene [String] the gene ID for KEGG data base
        @kegg_ID_pathway = {}      # kegg_IDs_pathway [Hash] the KEGG IDs and their corresponding description of the pathway
        self.annotate_GO
        self.annotate_kegg
    end

    # @!method no_params_method
    # Search for KEGG (IDs and corresponding descriptions of pathways) annotations for a give member.
    #   It needs the kegg gene ID to make the search
    # @return [void]
    def annotate_kegg
        return if @kegg_gene.empty?

        result = togo_search("kegg-genes", @kegg_gene)
        if result[0].key?("pathways") && !result[0]["pathways"].nil? && !result[0]["pathways"].empty? 
            @kegg_ID_pathway = result[0]["pathways"]
        end
        #puts "EL hash es: #{@kegg_ID_pathway}"
    end

    # @!method no_params_method
    # Search for GO (IDs and corresponding terms) annotations for a give member.
    #   It needs searchs only those GO terms referred to biological process.
    #   It also searchs for the kegg gene ID, which is needed later for kegg annotation.
    # @return [void]
    def annotate_GO
        result = togo_search("ebi-uniprot", self.uniprot_id, "/dr")
        if result[0].key?("GO") && !result[0]["GO"].nil? && !result[0]["GO"].empty?
            list_of_GOs = result[0]["GO"]
            list_of_GOs.each do |go|
                biological_process = go[1].match(/^[P]:/) #molecular functions (F), biological processes (P), and cellular components (C).
                if biological_process
                    id = go[0]
                    term = biological_process.post_match
                    @go_IDs_terms[id] = term
                else
                    next
                end
            end
        end

        if result[0].key?("KEGG") && !result[0]["KEGG"][0][0].nil? && !result[0]["KEGG"][0][0].empty?
            @kegg_gene = result[0]["KEGG"][0][0]
        end
        
    end

end