require 'json'

# @direct_interactors con attr_accesor y @network con get y set. Aunque @direct_interactors está en initialize
#PROBLEMA CON EL @network (se va a ir sobrescribiendo a medida que encuentre una nueva)
class Members

    attr_accessor :uniprot_id, :times_searched, :direct_interactors

    @@coexpresed_members= []
    @@all_members = Hash.new

    def initialize( uniprot_id: "IDXXX" )
        @uniprot_id = uniprot_id
        @times_searched = 0
        @direct_interactors = []
        @@all_members[@uniprot_id] = self
    end

    def gene_id=(gene_name) #SET gene_id
        @gene_id = gene_name
    end

    def gene_id #GET gene_id
        @gene_id
    end

    def set_network=(network)
        @network = network
    end
    
    def get_network
        @network
    end
    
    def self.read_from_file(filename) # Leer archivo y crear Members con cada ATG...
        coexpressed_file = File.open(filename, 'r')
        coexpressed_file.readlines.each do |line|
            locus_name=line.chomp
            address = "http://togows.dbcls.jp/entry/uniprot/#{locus_name}/accessions.json"
            response = RestClient::Request.execute(  #  or you can use the 'fetch' function we created last class
                method: :get,
                url: address)  
            result = JSON.parse(response.body)
            if result.is_a?(Array) && result.any?
                uniprot_id = result.first.first   # get the first uniprot id from list, there might be more than one
            else
                puts "No UniProt entry found for locus #{locus_name}. Please remove this entry from gene list"
                next
            end
            member = self.new(uniprot_id: uniprot_id)
            member.gene_id=(locus_name)
            @@coexpresed_members << member
        end
    end

    def find_interactors #Coge una instancia Member, busca sus interactores, y los crea (si no existen) añadiéndolos a @direct_interactors
        intact_address = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/'
        species = 'species:arabidopsis'
        formato = 'tab25'
      
        address = "#{intact_address}search/interactor/#{@uniprot_id}/?query=#{species}&format=#{formato}"
        puts "El miembro #{@uniprot_id} con #{address}"
        response = RestClient::Request.execute(
          method: :get,
          url: address,
          headers: {'Accept' => 'application/json'})
      
        if response.empty? #Si la rpta es en blanco, deja el array vacío
          #puts "#{member.gene_id} es respuesta en blanco"
          @direct_interactors = []
        else
            response.body.each_line do |line|
                values = line.chomp.split("\t")
                interactor_1 = values[0].match(/uniprotkb:([OPQ][0-9][A-Z0-9]{3}[0-9])/)
                interactor_2 = values[1].match(/uniprotkb:([OPQ][0-9][A-Z0-9]{3}[0-9])/)

                #Si mi interactor está en el campo 1 del tab25
                if !interactor_1.nil? && !interactor_1[1].include?(self.uniprot_id)

                    if @@all_members.include?(interactor_1[1]) # Si ya existe
                        @direct_interactors << @@all_members[interactor_1[1]]
                    else
                        interactor = self.class.new(uniprot_id: interactor_1[1]) 
                        @direct_interactors << interactor
                    end # Si todavía no existe

                #Si mi interactor está en el campo 2 del tab25
                elsif !interactor_2.nil? && !interactor_2[1].include?(self.uniprot_id)

                    if @@all_members.key?(interactor_2[1]) # Si ya existe
                        @direct_interactors << @@all_members[interactor_2[1]]
                    else # Si todavia no existe
                        interactor = self.class.new(uniprot_id: interactor_2[1])
                        @direct_interactors << interactor
                    end

                else # Si campo 1 y campo 2 fuesen mi query: este else es importante por si fuese la query A y fuese A-A
                    next
                end
            end
        end
    end

    def add_to_network #Coge los @direct_interactors de la instancia Member y 1)Define su @network como el del objeto Member 2)Los añade al Networks del Member
        unless @direct_interactors.empty?
            @direct_interactors.each do |interactor|
                unless self.get_network.network_members.key?(interactor.uniprot_id) #Si ya existe en la red del AT..., pasa. CREO QUE PODRÍAS PREGUNTAR SI ESTÁ ESE OBJETO INTERACTOR DIRECTAMENTE
                    interactor.set_network=(self.get_network)
                    self.get_network.add_member(interactor)
                end
            end
        end
    end

    def register_search
        @times_searched += 1
    end

    def self.all_coexpresed_members
        return @@coexpresed_members
    end

    def self.all_members
        return @@all_members
    end

    def eql?(other) 
        self.uniprot_id == other.uniprot_id if other.is_a?(Networks) # just to make sure we are comparing objects of the same class
    end

    def hash    # this code generates a hash code based on the attribute values
        # it is important for the correct fucntioning of hash-based collections (like Ruby's Hash),
        # since we are storing our netmembers in a hash, we do this so when we look for duplicates, we do it by id
        @uniprot_id.hash
    end
end