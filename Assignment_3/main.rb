#----------------- LOAD EXTERNAL LIBRARIES AND MODULES-------------------------------------------------------------------
require 'bio'
require 'net/http'
require  './EmblEntry.rb'
require './EmblProcessor.rb'

#------------------------------ INPUT FILES NAMES AS ARGUMENTS FROM COMMAND LINE ----------------------------------------

# Check for one common error: not specifying the correct number of files needed for the program to run
if ARGV.length != 1
    abort "Incorrect number of files passed. Gene file name must be specified: input list of genes"
end
  
# Check for second common error: incorrect usage, files in incorrect order or wrong name passed.
if ARGV[0] == "ArabidopsisSubNetwork_GeneList.txt"
    gene_file = ARGV[0]
else 
    abort "Incorrect order of files passed. Usage: ruby main.rb ArabidopsisSubNetwork_GeneList.txt"
end


# ------------ PUBLIC INSTANCE METHODS ---------------------------------------

# Gets a response from a given URL
#
# @param uri_str [String] URL to make the get request to 
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

# Programmatic access to NCBI E-Utils,
#
# @param database [String] the database to retrieve gene's information
# @param file_format [String] type of format for the response
# @param gene_locus [String] the gene's locus name to retrieve info about.
def ncbi_fetch(database, file_format, id)
    url = "http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=#{database}&format=#{file_format}&id=#{id}"
    #puts url
    res = fetch(url)
    return res
end


# Read all genes (locus names in the file) and create a comma-separated string of them. Check for correct gene format and species
#
# @param filename [String] Name of the file with the locus names of the genes to retrieve sequence's info about
# return gene_ids [Array<String>] array with all gene ids
def read_from_file(filename)
    gene_list = []
    gene_file = File.open(filename, 'r')
    gene_file.readlines.each do |line|  
        locus_name=line.chomp
        if locus_name !~ /A[Tt]\d[Gg]\d\d\d\d\d/ # Checking proper arabidopsis thaliana gene locus
            abort "Locus name #{locus_name} does not meet the correct format. Please define locus names as AT0g00000"
        end
        gene_list << locus_name 
    end
    gene_file.close
    gene_ids = gene_list.join(",")  # All gene IDs will be retrieved at once from embl
    return gene_ids
end

# Given a string with all cordinates it splits and finds start and end coordinates
#
# @param coordinates_string [String] string of coordinates to retrieve start and end coordinates
# return [Array<Integer>] retrieved coordinates as a two-element array
def get_coordinates(coordinates_string)
    position = coordinates_string.match(/\d+\.\.\d+/)[0]  # To match coordinates
    start_coordinate = position.split('..')[0].to_i
    end_coordinate = position.split('..')[1].to_i
    return[start_coordinate, end_coordinate]
end

# Given an specific motif (sequence pattern) get its start and end coordinates of any match in a given sequence
#
# @param sequence [String] sequence string to find the specific pattern
# @param motif [String] string pattern 
# return matches [Array<String>] array with all matches (found motifs) and its start and end coordinates (referenced to the sequence string param)
def find_motifs (sequence, motif)
    matches = []
    position = 0    # Start position to search
    while (match = sequence.match(motif, position))
        match_sequence = match[0]  
        start_motif = match.begin(0)    # Position of the first character to match
        end_motif = match.end(0)        # Position of the second character to match
        matches << {   
          sequence: match_sequence,
          start: start_motif.to_i,
          end: end_motif.to_i
        }

        position = end_motif - 5    # Move the searching position for next match, if existing
    end
    return matches
end

# Write a GFF3 file
#
# @param gff [Bio::GFF::GFF3] Bio::GFF object to be printed out into a file
# @param coordinates_to_use [String] Reference system to calculate motif's coordinates respect with 

def write_gff(gff, coordinates_to_use)
    if coordinates_to_use == 'chromosomal coordinates'
        file = 'AT_repeats_chromosomal.gff'
    else
        file = 'AT_repeats_sequence.gff'
    end
    File.open(file, 'w') do |file|
      file.puts gff.to_s
    end
end
  
#-------------------------------------------MAIN CODE ------------------------------------------------------------------

# Task 1

puts "Processing #{gene_file} file, this might take a while..."

gene_ids = read_from_file(gene_file) # Array of gene ids
response_body = ncbi_fetch(database = 'ensemblgenomesgene',file_format = 'embl', id = gene_ids) # One NCBI search for all genes 

filename = 'AT_sequences.embl' # Name of the embl file

output_file = File.open(filename, 'w')
output_file.write(response_body)
output_file.close


search_positive = Bio::Sequence::NA.new("cttctt")   # Motif to look in positive strand
search_complementary = Bio::Sequence::NA.new("aagaag")  # Motif to look for in negative strand
REGEX_POSITIVE = Regexp.new(search_positive.to_re)
REGEX_COMPLEMENT = Regexp.new(search_complementary.to_re)


# Tasks 2, 3, 4a and 4b

eobject_processor = EmblProcessor.new(filename)
gff_seq = eobject_processor.load_to_gff('sequences coordinates')    # Create the GFF3 file referred to the gene's sequence coordinates
write_gff(gff_seq, 'sequences coordinates')

# Tasks 5 and 6

eobject_processor_chr = EmblProcessor.new(filename)
gff_chr = eobject_processor_chr.load_to_gff('chromosomal coordinates')  # Create the GFF3 file referred to whole chromosome
write_gff(gff_chr,'chromosomal coordinates')



