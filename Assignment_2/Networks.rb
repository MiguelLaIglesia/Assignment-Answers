class Networks

  attr_accessor :network_members

  @@number_of_networks = 0
  @@all_networks = []

  def initialize
    @network_members = []
    @@all_networks << self
  end

  def add_member(net_member)
    @network_members << net_member #unless @network_members.include?(net_member), no hace falta, ya está en la función de recursividad
  end

  def self.get_all_networks
    return @@all_networks
  end

  def merge_networks
    @@all_networks.each do |network|
      next if network == self
      if network.network_members.intersection(@network_members) #Si hay miembros en comun entre ambas redes
        #puts "OJOOOOOOOOOOO"
        #puts "#{network.network_members} coincide con #{@network_members}"
        @network_members = @network_members | network.network_members # Se suman los miembros de ambas redes sin repetir, almacenándose en la red que llama a esta función

        network.network_members.each do |member| #Actualiza el valor de @network de los miembros absorbidos
          member.set_network=(self)
        end
        @@all_networks.delete(network) # Se elimina de @all_networks aquella que ha sido absorbida
      end
    end
  end

end