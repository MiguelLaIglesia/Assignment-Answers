# REVISAR NCBI LISTA IDs
# exon_sequence = entry_sequence.subseq(start_exon, end_exon): si el exon no cae en seq es que no está?
# COMO EVITAR OVERLAP DE MISMAS FEATURES
# REVISAR QUE TODAS LAS LAS FEATURES ANOTADAS TENGAN UN FORMATO ESTÁNDAR EN LA PAGINA ESA
# EL ID EN ATTRIBUTES ES COMO QUIERAS PERO IGUAL PARA TODOS Y CON UNA ENUMERACIÓN NO?ç
# Ver SI EL SOURCE ES ESO O EL "annotated by Araport11", yo creo que eso
# DOMINIOS CTTCTT QUE SE REFIEREN A LA MISMA SECUENCIA, PERO QUE ESTÁN EN 2 EXONES DISTINTOS, CÓMO CONSIDERARLO?
# UNA SOLA LINEA DEL GFF Y LUEGO EN ATRIBUTES LE METES TODOS LOS EXONES QUE LO TIENEN?
#Añadir una linea de metadata??
# hacer un filtro por "note" del exon que incluya el gen de la lista ? realmente no es necesario pues los de ccs de la secuencia son los q se miraran y por
#tanto son del gen de la lista.


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

# Gets coordinates from string
def get_coordinates(coordinates_string)
    position = coordinates_string.match(/\d+\.\.\d+/)[0]

    start_coordinate = position.split('..')[0].to_i
    end_coordinate = position.split('..')[1].to_i
    return start_coordinate, end_coordinate
end


def find_motifs (sequence, motif)
    matches = []
    position = 0
  
    while (match = sequence.match(motif, position))
        match_sequence = match[0]
        start_motif = match.begin(0)
        end_motif = match.end(0)
  
        matches << {
          sequence: match_sequence,
          start: start_motif.to_i,
          end: end_motif.to_i
        }
  
        # Move the position forward, "- 5" allows to detect a domain "CTTCTTCTT", and that will be detected as 2
        position = end_motif - 5
    end
  
    return matches
end



#-------------------------------------------MAIN CODE ------------------------------------------------------------------
                
                # DESACTIVADO TODO, ACTIVAR PARA ENTREGA
                # 1:  Using BioRuby, examine the sequences of the ~167 Arabidopsis genes from the last assignment by retrieving them from whatever database you wish #
                # ArabidopsisSubNetwork file -> each gene locus (lines) -> search for embl file -> add to 'AT_sequences.embl'

                #puts "Processing #{gene_file} file, this might take a while..."

                # Create string with gene IDs from file
                #gene_ids = read_from_file(gene_file)

                # Only one search with list of IDs
                #response_body = ncbi_fetch(database = 'ensemblgenomesgene',file_format = 'embl', id = gene_ids)

                #output_file = File.open('AT_sequences.embl', 'w')
                #output_file.write(response_body)
                #output_file.close

                # Comprobado que da los 168

# 2: Loop over every exon feature, and scan it for the CTTCTT sequence #

embl_file = Bio::FlatFile.auto('AT_sequences.embl')

search_positive = Bio::Sequence::NA.new("cttctt")
search_complementary = Bio::Sequence::NA.new("aagaag")

regex_positive = Regexp.new(search_positive.to_re)
regex_complementary = Regexp.new(search_complementary.to_re)

i = 0


source = ""
entries = []
k=0
h=0

# Loop each entry in EMBL file
embl_file.each_entry do |entry|

    #break if i > 20
    #i += 1

    next unless entry.accession # Lack of accesion is suspicious

    #puts
    #puts "NEW ENTRY #{entry.definition}"
    
    # Extract sequence chromosome coordinates from entry.definition
    start_entry, end_entry = get_coordinates(entry.definition)

    # Anotate sequence chromosome coordinates
    entry_sequence = entry.to_biosequence

    f1 = Bio::Feature.new('chromosomal_coordinates', entry.definition) # 'feature type' , 'position'
    f1.append(Bio::Feature::Qualifier.new('start', start_entry))
    f1.append(Bio::Feature::Qualifier.new('end', end_entry))
    entry_sequence.features << f1

    # LOOP FEATURES
    entry.features.each do |feature|
        featuretype = feature.feature  
        position = feature.position 
        qual = feature.assoc

        if featuretype == 'source'
            source = qual['db_xref']
        end

        # To annotate Arabidopsis thaliana gene as a new feature. Tests if it was not annotated before and if the format of coordinates are correct.
        if featuretype == 'gene' && feature.position.match(/^(complement\()?(\d+\.\.\d+)(\))?$/)

            # Test if that feature already exists
            feature_exists = entry_sequence.features.any? { |feature| feature.feature == 'gene_from_list'}
            next if feature_exists

            # New feature with AT gene name ("AT5G48300")
            f1 = Bio::Feature.new('gene_from_list', qual['gene']) # 'feature type' , 'position'
            entry_sequence.features << f1
        end
        
        # Tests if feature is "exon" and if the coordinates do not correspond to another fragment
        unless featuretype == 'exon' && position.match?(/^(complement\()?(\d+\.\.\d+)(\))?$/)
            next
        end

        # Get exon ID
        exon_id = qual["note"]

        # Get exon coordinates
        start_exon, end_exon = get_coordinates(position)
        
        # Get exon sequence from its coordinates, and test if that exon is inside the entry sequence to continue or not
        exon_sequence = entry_sequence.subseq(start_exon, end_exon)
        next if exon_sequence.nil? || exon_sequence.empty?

        # Decides wether to use positive or complement regular expression
        if position.include?('complement')
            regex = regex_complementary
            strand = '-'
            coordinates_format = "complement(%s)" # %s indicates the substituient, in this case, coordinates range
        else
            regex = regex_positive
            strand = '+'
            coordinates_format = "%s"
        end

        # Find motifs in exon sequence and retrieves the positions (BUT IN EXON SEQUENCE). It works for several matches/motifs in the exon.
        matches_motifs = find_motifs(exon_sequence, regex)
        
        next if matches_motifs.empty? || matches_motifs.nil?
        
        # For each matched motif in exon sequence
        matches_motifs.each do |match|

            # First and last position of motif, BUT IN EXON SEQUENCE
            start_match = match[:start]
            end_match = match[:end]

            # Conversion from exon positions to entry sequence positions
            start_motif = start_exon + start_match # +1 por ccs del match, -1 por la suma de exon y match
            end_motif = start_exon + end_match - 1 # -1 por la suma de exon y match
            #puts "Motif #{entry_sequence.subseq(start_motif, end_motif)} found at position #{start_motif} to #{end_motif} in #{exon_id}"
           
            coordinates = coordinates_format % "#{start_motif}..#{end_motif}"

            # Tests if that feature has already been annotated, as several exons have same ctttcttt positions
            feature_exists = entry_sequence.features.any? { |feature| feature.feature == 'cttctt_repeat' && feature.position == coordinates }
            next if feature_exists

            k = k + 1

            # Anotates new feature
            f1 = Bio::Feature.new('cttctt_repeat', coordinates ) # 'feature type' , 'position'
            f1.append(Bio::Feature::Qualifier.new('start', start_motif))
            f1.append(Bio::Feature::Qualifier.new('end', end_motif))
            f1.append(Bio::Feature::Qualifier.new('strand', strand))
            #f1.append(Bio::Feature::Qualifier.new('score', '.')) # In this case, there is no score value
            #f1.append(Bio::Feature::Qualifier.new('phase', '.')) # This is only set when feature type is "CDS", not "exon"
            f1.append(Bio::Feature::Qualifier.new('source', source))
            f1.append(Bio::Feature::Qualifier.new('SO_Name', 'repeat_region'))
            entry_sequence.features << f1


        end

    end

    # Keeping all entry objects in a list for later
    entries << entry_sequence

end

puts
puts

# GFF FILE
gff = Bio::GFF::GFF3.new

count = 1 # To annotate attributes

entries.each do |entry|

    # To test if that entry has a repetitive motif annotated, this value True/False will be used later
    feature_exists = entry.features.any? { |feature| feature.feature == 'cttctt_repeat' }

    entry.features.each do |feature|

        featuretype = feature.feature
        qual = feature.assoc
        
        # For those without repetitive motif, puts it gene locus ("AT5G48300")
        if featuretype == 'gene_from_list' && !feature_exists
            puts "There are no repetitive cttctt regions found for #{feature.position}"
            break
        end

        next unless featuretype == 'cttctt_repeat'
        h = h + 1
        
        #puts "Chr#{entry.entry_id} #{qual['source'].class} #{qual['SO_Name'].class} #{qual['start'].class} #{qual['end'].class} #{qual['strand'].class} ID=repeat_region_#{count};Name=CTTCTT_motif"
        
        attributes = [{"ID" => "repeat_region_#{count}", "Name" => "CTTCTT_motif"}]
        
        gff.records << Bio::GFF::GFF3::Record.new(
            "Chr#{entry.entry_id}",      # seqID
            qual['source'],     # source
            qual['SO_Name'],    # feature type
            qual['start'],            # start
            qual['end'],          # end
            nil,          # score
            qual['strand'],          # strand
            nil,          # phase
            attributes[0] # attributes
         )
        count += 1
    end
end



puts
puts k
puts h

#gff.records.each do |feature|
#    puts feature.class
#end

File.open('output_1.gff', 'w') do |file|
    file.puts gff.to_s
end