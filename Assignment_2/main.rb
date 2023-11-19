require 'rest-client'
require 'json'
require './Members'
require './Networks'

Members.read_from_file('ArabidopsisSubNetwork_GeneList.txt')

#Comprobado, el AT... coincide con el uniprot_id
#Members.all_coexpresed_members.each do |member|
#  puts 'New'
#  puts member.gene_id
#  puts member.uniprot_id
#end

#puts Members.all_members


def recursive_search(member, network, depth)

  #return 1 if depth > 2 PARA CORTAR LA RECURSIVIDAD: ya será 3 cuando hagas la nueva busqueda en 2
  if depth == 1
    puts "Estamos en nivel #{depth}"
    member.find_interactors if member.direct_interactors.empty? # Buscalos solo si antes no los buscaste o si no es una string de "Not available"
    member.add_to_network
    if member.direct_interactors.empty? || member.direct_interactors.is_a?(String)
      puts "Interacciones para #{member.uniprot_id} no encontradas en IntAct, corte de la recursividad"
      return 1
    end
    return recursive_search(member.direct_interactors, network, depth = depth + 1)

  elsif depth == 2
    puts "Estamos en nivel #{depth}"
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

  else
    puts "Hemos pasado el depth 2"
    return 1
  end


  #puts "Los direct_interactors son #{member.direct_interactors}"
  #puts "Los miembros de la red son #{network.network_members}"
  #puts
  #puts

end

Members.all_coexpresed_members.each do |member|
  puts
  puts
  puts "ANALIZANDO NUEVO MIEMBRO con #{member.gene_id} y #{member.uniprot_id}"
  network = Networks.new()
  member.set_network=(network)
  network.add_member(member)
  recursive_search(member, network, depth=1)
  puts "Este miembro tiene la red #{network} con miembros:"
  network.network_members.each do |_key, value|
    puts value.uniprot_id
  end
    
  #puts "Para #{member.uniprot_id} con #{member.gene_id} tenemos la red #{network} con miembros #{network.network_members}"
  #puts
  #puts
end