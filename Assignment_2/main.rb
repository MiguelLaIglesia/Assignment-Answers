#----------------- REQUIRED MODULES OR CLASSES ------------------#

require 'rest-client'
require 'json'
require './Members'
require './Networks'
require './AnnotatedMembers'

#----------------- MAIN FUNCTIONS ------------------#

# Request programatic access to URL
# @param address [String] the URL to search
# @return [void]

def rest_api_request(address)
  RestClient::Request.execute(
    method: :get,
    url: address,
    headers: {'Accept' => 'application/json'})
end

# Extracts relevant information from tab25 format, i.e., the detection method used for that interaction
# @param element [String] the field of tab25 file
# @return [String] the information extracted

def extract_xref(element)
  match_data = element.match(/:(\S+)(?:\(|$)/)
  val = match_data[1].gsub(/\A"/, '').gsub(/"\z/, '') if match_data
  return(val)
end

# Search in TOGO REST API to access databases' information retrieving a JSON output file
# @param database [String] the database required
# @param query [String] the query (i.e. protein identifier or kegg pathway) to search
# @param field [String] the directory to index in the URL search
# @return [String] the parsed JSON ouput file with the required information

def togo_search(database,query,field="")
  togo_address = "http://togows.dbcls.jp/entry/#{database}/#{query}#{field}.json"
  togo_response = rest_api_request(togo_address)
  result = JSON.parse(togo_response.body)
  return result
end

#----------------- INPUT FILES IN COMMAND LINE ------------------#

# Testing number of arguments equal to 2
if ARGV.length != 2
  abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end

# Testing proper name and order of arguments
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt" && ARGV[1] == "Final_report.txt"
  input_gene_list = ARGV[0]
  output_report_file = ARGV[1]
else 
  abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt Final_report.txt"
end


#----------------- MAIN CODE ------------------#

puts "Processing #{input_gene_list} file, this might take a while..."

# Calling function for creating members from each gene (line) in input file
LuMike_objects::Members.read_from_file(input_gene_list)

# Parameters set for this assignment
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
SPECIES = 'species:arabidopsis' 
TAB25_FORMAT = 'tab25'
DEPTH = 2

# Each member from the coexpressed genes list is iterated
LuMike_objects::Members.all_coexpresed_members.each do |gene|
  puts
  puts "Analyzing new gene #{gene.gene_id} with Uniprot ID #{gene.uniprot_id}"
  network = LuMike_objects::Networks.new # Seed network
  network.add_member(gene) # Add gene to network
  network.recursive_search([gene], DEPTH) # Recursive search for the gene
  network.merge_with_common # Merge that network to the ones already existing
  puts "This member is part of network #{network} with #{network.network_members.length} miembros"
end

# Calling method to delete networks with only one member
LuMike_objects::Networks.reduce_networks

# Each network is annotated with KEGG and GO
LuMike_objects::Networks.all_networks.each do |network|
  network.annotate_network
end

# Report the results to a file
File.open(output_report_file, 'w') do |file|

  file.puts "This code was created by Miguel La Iglesia Mirones and Lucía Muñoz Gil"
  file.puts "Bioinformatics Programming Challenges Course at MSc Computational Biology (UPM)"
  file.puts "November 2023"
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "Protein - protein interaction networks (interactome) final report for genes detailed" 
  file.puts "at #{input_gene_list} file"
  file.puts
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "GLOBAL REPORT FOR ALL PROTEIN-PROTEIN INTERACTION NETWORKS"
  file.puts
  file.puts "Total number of nodes: #{LuMike_objects::Members.number_of_members}"
  file.puts "Number of nets: #{LuMike_objects::Networks.get_number_of_nets}"
  file.puts "Genes from list not included in any network: #{LuMike_objects::Networks.nodes_without_net}"
  file.puts
  file.puts "---------------------------------------------------------------------------------------"
  file.puts
  file.puts "FEATURES OF EVERY INTERACTION NETWORK"
  LuMike_objects::Networks.all_networks.each_with_index do |network, idx|
    file.puts
    file.puts "Network ID #{idx + 1} (#{network}):"
    file.puts "  Number of members: #{network.network_members.length}"
    file.puts "  Genes from file part of this network:"
    
    network.network_members.each do |member|
      if member.gene_id
        file.puts "    Gene ID: #{member.gene_id}, Uniprot ID for coded protein: #{member.uniprot_id}" 
      end
    end

    file.puts
    file.puts "  KEGG annotations in network #{idx + 1}:"
    network.network_KEGGs.each do |key_, value|
      file.puts "    KEGG ID: #{key_}   Pathway: #{value}"
    end

    file.puts
    file.puts "  GO annotations in network #{idx + 1}:"
    network.network_GOs.each do |key_, value|
      file.puts "    GO ID: #{key_}   Term: #{value}"
    end

  end 
end