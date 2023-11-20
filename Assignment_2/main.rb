require 'rest-client'
require 'json'
require './Members'
require './Networks'

#------------------MAIN FUNCTIONS --------------------------------------------

def rest_api_request(adress)
  RestClient::Request.execute(
    method: :get,
    url: adress,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute"
end

def extract_xref(element) # syntax of TAB25 <XREF><VALUE>(<DESCRIPTION>)
  # match regex to get <VALUE>
  match_data = element.match(/:(\S+)(?:\(|$)/)  
  # get the actual value
  val = match_data[1] if match_data
  return(val)
end

#-----------------------------------------------------------------------------

Members.read_from_file('ArabidopsisSubNetwork_GeneList.txt')


# Parameters for this assignment
INTACT_BASE_ADDRESS = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
SPECIES = 'species:arabidopsis' 
TAB25_FORMAT = 'tab25'
DEPTH = 2




Members.all_coexpresed_members.each do |gene|
  puts
  puts
  puts "ANALIZANDO NUEVO MIEMBRO con #{gene.gene_id} y #{gene.uniprot_id}"
  
  network = Networks.new  # create the new network
  #gene.set_network=(network)
  network.add_member(gene)
  network.recursive_search([gene], DEPTH) # search and assign all the interactors found to this net

  puts network.create_and_merge

  puts "Este miembro tiene la red #{network} con #{network.network_members.length} miembros"

  #network.network_members.each do |element|
  #  puts element.uniprot_ids
  #end
end


File.open('Resultado.txt', 'w') do |file|

  file.puts "------------FIN--------------"
  file.puts "TODAS LAS REDES CREADAS SON"
  Networks.all_networks.each do |network|
    file.puts
    file.puts "Red #{network} con miembros:"
  
    network.network_members.each do |member|
      if member.gene_id
        file.puts "#{member.uniprot_id}, aqui está el gen: #{member.gene_id}" 
      else
        file.puts member.uniprot_id
      end
    end
  end 
  file.puts
  file.puts "Numero de miembros en total: #{Members.all_members.length}"
  file.puts
  file.puts "Número de redes en total: #{Networks.all_networks.length}"
end