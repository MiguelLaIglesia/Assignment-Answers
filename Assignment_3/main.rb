#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------

require 'bio'
require 'net/http'

#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM  COMMAND LINE ---------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 1
    abort "Incorrect number of files passed. Two files names must be specified: input list of genes and name for final_report"
end
  
# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt"
    gene_file = ARGV[0]
else 
    abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt"
end





#------------------------ PUBLIC INSTANCE METHODS -----------------------------------------------------------------------



# Get response from URL
def fetch(uri_str) 
    address = URI(uri_str)  
    response = Net::HTTP.get_response(address)
    case response
      when Net::HTTPSuccess then
        return response.body
      else
        raise Exception, "Something went wrong... the call to #{uri_str} failed; type #{response.class}"
        response = false
        return response
    end
end

# Get response with NCBI E-Utils, given a database, file format and gene locus
def ncbi_fetch(database, file_format, id)
    url = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=#{database}&format=#{file_format}&id=#{id}"
    #puts url
    res = fetch(url)
    return res
end


# Transform locus names file into string of coma separated locus names, while checking A. thaliana format
def read_from_file(filename)
    gene_list = []
    gene_file = File.open(filename, 'r')
    gene_file.readlines.each do |line|
        
        locus_name=line.chomp

        if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Check for correct format of gene locus name
            abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
        end

        gene_list << locus_name
        
    end

    gene_file.close

    gene_ids = gene_list.join(",")

    return gene_ids
end

#-------------------------------------------MAIN CODE ------------------------------------------------------------------

# 1:  Using BioRuby, examine the sequences of the ~167 Arabidopsis genes from the last assignment by retrieving them from whatever database you wish
# ArabidopsisSubNetwork file -> each gene locus (lines) -> search for embl file -> add to 'AT_sequences.embl'
# DESACTIVADO TODO, ACTIVAR PARA ENTREGA

#puts "Processing #{gene_file} file, this might take a while..."

# Create string with gene IDs from file
#gene_ids = read_from_file(gene_file)

# Only one search with list of IDs
#response_body = ncbi_fetch(database = 'ensemblgenomesgene',file_format = 'embl', id = gene_ids)

#output_file = File.open('AT_sequences.embl', 'w')
#output_file.write(response_body)
#output_file.close

# Comprobado que da los 168

# 2: Loop over every exon feature, and scan it for the CTTCTT sequence

embl_file = Bio::FlatFile.auto('AT_sequences.embl')
#puts embl_file.class  # Bio::FlatFile

#entry = embl_file.next_entry

search = Bio::Sequence::NA.new("cttctt")
regex = Regexp.new(search.to_re)

i=0

embl_file.each_entry do |entry|

    break if i > 2
    i += 1

    next unless entry.accession 

    puts
    puts "NUEVA ENTRY #{entry.definition} #{entry.class}"
    # Extract sequence coordinates from entry.definition
    entry_definition = entry.definition
    match = entry_definition.match(/\d+\.\.\d+/)
    if match
        coordinates = match[0]
        start_entry = coordinates.split('..')[0]
        end_entry = coordinates.split('..')[1]
    end

    entry_sequence = entry.to_biosequence
    f1 = Bio::Feature.new('seq_coordinates',coordinates) # feature type and position
    f1.append(Bio::Feature::Qualifier.new('start', start_entry))
    f1.append(Bio::Feature::Qualifier.new('end', end_entry))
    entry_sequence.features << f1

    #puts entry_sequence.length
    #entry_sequence = entry.naseq
    #entry_sequence = entry.sequence

    #puts entry_sequence.class
    #entry_id = entry_sequence.entry_id
    #entry_end = entry_sequence.end

    #puts "#{entry_id}"
    #entry_position = entry.position
    #puts  entry_position
    #puts sequence
    #puts sequence.seq
    #seq = Bio::Sequence::NA.new("atgcatgcaaaa")

    #SACAR COORD TOTALES DE ENTRY DEFINITION Y DE AHI BUSCAR EXON EN LA SECUENCIA (CON SU POSITION) Y VER SI HAY CTTCTT
    entry.features.each do |feature|
        
        featuretype = feature.feature

        next unless featuretype == 'exon'

        position = feature.position # F23J3_4:28298..28420
        puts position

        # DIFERENCIAR ENTRE COMPLEMENT O NO PARA HCER SPLIT DISTINTO Y TRABAJAR CON COORDENADAS INTERNAS O TOTALES
        # CREAR OBJETO SECUENCIA DEL EXON CON START Y END, YO DIRIA YA REFERIDA A TOTALES (ENTONCES LA + TRANSFORMAR A TOTALES, EASY)
        # Y DE AHI BUSCAR EN LA SECUENCIA SI HAY CTTCTT CON EL MATCH

        #start_exon = position.split('..')[0]
        #end_exon = position.split('..')[1]

        #puts start_exon
        #puts end_exon
        

        #puts featuretype
        #if featuretype == 'seq_coordinates' # estan conectados entry y entry_sequence
        #    
        #    puts "AAAAAAAAAAA #{feature.position}"
        #end
        
        # OBTENER COORDENADAS INTERNAS DE LA SECUENCIA
        #if feature.feature == "source"
        #    puts feature 
        #    puts feature.position
        #end

        #qual = feature.assoc
        #puts "#{feature} #{feature.position} #{qual}"
        
        #exon_seq = seq.subseq(3,8)

        #match = exon_seq.seq.match(regex)

        #qual = feature.assoc            # feature.assoc gives you a hash of Bio::Feature::Qualifier objects 
        #puts qual
        #puts position
        #puts "NUEVA FEATURE #{qual}"                                # i.e. qualifier['key'] = value  for example qualifier['gene'] = "CYP450")
        #next unless qual['exon']    # this is an indication that the feature is a transcript

        #exon = qual["note"]
        #puts "#{exon}"
        #puts exon
        # collects gene name and exon positions and joins it into a string
        #gene_info = [qual['gene'], qual['exon']].compact.join(', ')
        #puts "aaaaa #{exon}"
    end
end









# LOOP FILE ENTRIES (MÁS EN JUPYTER)
#embl_file.each_entry do |entry|
#    next unless entry.accession 
#    puts entry.class
#    puts "The record is #{entry.definition} "
#  end
