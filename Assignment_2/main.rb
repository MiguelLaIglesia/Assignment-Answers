require 'rest-client'
require 'json'
require './Members'
require './Networks'

#CAMBIOS: funcion recursiva entera, bucle genes y file entero

#------------------MAIN FUNCTIONS --------------------------------------------

def rest_api_request(adress)
  RestClient::Request.execute(
    method: :get,
    url: adress,
    headers: {'Accept' => 'application/json'})  # use the RestClient::Request object's method "execute",
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


def recursive_search(member, network, depth)

  return 1 if depth > 2

  if depth == 1
    #puts "Estamos en nivel #{depth}"
    member.find_interactors if member.direct_interactors.empty? # Buscalos solo si antes no los buscaste o si no es una string de "Not available"
    member.add_to_network
    if member.direct_interactors.empty? || member.direct_interactors.is_a?(String)
      #puts "Interacciones para #{member.uniprot_id} no encontradas en IntAct, corte de la recursividad"
      return 1
    end
    return recursive_search(member.direct_interactors, network, depth = depth + 1)

  else
    #
    #puts "Estamos en nivel #{depth}"
    list_of_interactors = []
    #puts "longuitud #{member.length}"
    member.each do |one_member|
      #puts "Este es #{one_member} con #{one_member.uniprot_id}"
      one_member.find_interactors if one_member.direct_interactors.empty?
      one_member.add_to_network
      list_of_interactors += one_member.direct_interactors
    end
    #list_of_interactors.each do |list| #Funciona
    #  puts list.uniprot_id
    #end
    return recursive_search(list_of_interactors, network, depth = depth + 1)
  end


  #puts "Los direct_interactors son #{member.direct_interactors}"
  #puts "Los miembros de la red son #{network.network_members}"
  #puts
  #puts

end

Members.all_coexpresed_members.each do |gene|
  #puts
  #puts
  #puts "ANALIZANDO NUEVO MIEMBRO con #{gene.gene_id} y #{gene.uniprot_id}"

  network = Networks.new()
  gene.set_network=(network)
  network.add_member(gene)

  recursive_search(gene, network, depth=1)

  gene.get_network.merge_networks


  #puts "Este miembro tiene la red #{network} con miembros:"

  #network.network_members.each do | net_member|
  #  puts net_member.uniprot_id
  #end

end

File.open('Resultado.txt', 'w') do |file|
  file.puts
  file.puts
  file.puts "------------FIN--------------"
  file.puts "TODAS LAS REDES CREADAS SON"
  Networks.get_all_networks.each do |network|
    puts
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
  file.puts "Número de redes en total: #{Networks.get_all_networks.length}"
end