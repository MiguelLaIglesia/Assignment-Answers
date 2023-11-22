# == Networks
#
# This is a representation of protein interaction networks
# with the network members, network GO and KEGG annotations
#
# == Summary
# 
# This can be used to represent different interactors that belong to the same network
#


class Networks

  # Get/Set the members of the network
  # @!attribute [rw]
  # @return [array<Members>] network members

  attr_accessor :network_members

  # Get/Set the number of networks
  # @!attribute [rw]
  # @return [Integers] total number of networks

  attr_accessor :number_of_networks

  # Get/Set the network GO annotations
  # @!attribute [rw]
  # @return [array] network GO from all network members

  attr_accessor :network_GOs

  # Get/Set the network KEGG annotations
  # @!attribute [rw]
  # @return [array] network KEGG from all network members

  attr_accessor :network_KEGGs

  @@total_networks = 0
  @@nodes_without_net = 0
  @@all_networks = []
  
  # Create a new instance of Networks
  # @return [Networks] an instance of Networks

  def initialize
    @network_members = []       # network_memners [array<Members>] the members belonging to the network
    @@total_networks += 1       # total_networks [Integer] the total number of networks
    @network_GOs = {}           # network_GOs [Hash] the GO IDs and terms annotated for that network
    @network_KEGGs = {}         # network_KEGGs [Hash] the KEGG IDs and terms annotated for that network
    @@all_networks << self      # all_networks [array<Networks>] the list containing all networks
  end

  # @!method no_params_method
  # Get all networks
  # @return [array<Networks>] the total of networks

  def self.all_networks
    return @@all_networks
  end

  # @!method no_params_method
  # Get total number of networks
  # @return [Integer] the total number of networks

  def self.get_number_of_nets
    @@total_networks
  end

  # @!method no_params_method
  # Get all nodes without a network
  # @return [array<Members>] the members without a network

  def self.nodes_without_net
    @@nodes_without_net
  end

  # Adds a new member for an existing network
  # @param net_member [<Members>] the member to add
  # @return [void]

  def add_member(net_member)
    @network_members << net_member
    #@network_GOs.merge!(net_member.go_IDs_terms)
    #@network_KEGGs.merge!(net_member.kegg_ID_pathway)
  end

  # @!method no_params_method
  # Tests if the number of members for a network is lower than 2
  # @return [void]

  def only_one_member?
    @network_members.length < 2
  end

  # @!method no_params_method
  # Deletes those networks with only one member
  # @return [void]

  def self.reduce_networks
    before_reduction = @@all_networks.length
    @@all_networks = @@all_networks.select { |network| !network.only_one_member? }

    @@nodes_without_net = before_reduction - @@all_networks.length
    @@total_networks -= @@nodes_without_net
  end

  # @!method no_params_method
  # Annotates a network considering GO and KEGG annotations of its members
  #    It merges members these annotations so that no duplicates are stored
  # @return [void]

  def annotate_network

    @network_members.each do |member| 

      if !member.go_IDs_terms.nil? && !member.go_IDs_terms.empty?
        member.go_IDs_terms.each do |key, value|
          @network_GOs[key] = value unless @network_GOs.key?(key) # Adding GO if it is not already stored
        end
      end

      if !member.kegg_ID_pathway.nil? && !member.kegg_ID_pathway.empty?
        member.kegg_ID_pathway.each do |key, value|
          @network_KEGGs[key] = value unless @network_KEGGs.key?(key) # Adding KEGG if it is not already stored
        end
      end

    end
  end


  # Adds direct interactors of a member to a network
  # @param net_member [<Members>] the member to add its direct interactors to network
  # @return [void]

  def add_interactors_to_network(net_member)
    net_member.direct_interactors.each do |interactor|
      unless @network_members.include?(interactor) # Adding member if it is not already in network
            #interactor.set_network=(self)
            add_member(interactor)
      end
    end 
  end


  # Recursive search for interactors in a network
  # @param found_proteins [Array<Members] the list of members to search for its interactors in the network
  # @return [void]

  def recursive_search(found_proteins, depth)
    return if depth <= 0
  
    list_of_interactors = []
    found_proteins.each do |protein|
      protein.find_interactors if protein.direct_interactors.empty?
      next if protein.direct_interactors.empty? || protein.direct_interactors.is_a?(String)
  
      add_interactors_to_network(protein)
      list_of_interactors += protein.direct_interactors # Interactors are stored for calling recursivity later
    end
    
    recursive_search(list_of_interactors, depth - 1)
  end

  # Merges two networks that have common members between them
  #    It merges members that belong to each of both networks. One of the networks is deleted, the other is the one merged
  # @param other_network [<Networks>] the network to merge with the one that is calling the method
  # @return [void]

  def merge_with(other_network)

    common_members = @network_members & other_network.network_members # Tests if there are members in common
    if common_members.any?
      @network_members |= other_network.network_members
      @@all_networks.delete(other_network)  # Deleting old net: two were merged in one
      @@total_networks -= 1
    end
  end

  # @!method no_params_method
  # Merges networks with common members considering all existing networks.
  #   It calls merge_with to merge two networks with common interactors.
  # @return [void]

  def merge_with_common

    nets_with_common_members = @@all_networks.select do |existing_net|
      existing_net.network_members.any? { |member| self.network_members.include?(member)}
    end

    if nets_with_common_members.size > 1
      nets_with_common_members.each do |common_net|
        self.merge_with(common_net) unless common_net == self
      end
    end

  end

end