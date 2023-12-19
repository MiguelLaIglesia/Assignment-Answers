require 'bio'

#PROBLEMAS
# PARAMETROS BLAST
# SOLUCIONAR GUIONES
# HIT ID?

# opcion traducir o no desde linea de comandos o como parametro del script
# que decida que tipo de blast en función del tipo de secuencia

# -------------- MAIN FUNCTIONS -------------------------------------------------

# Check content of the fasta files of the proteomes, the sequences are genomic sequences or protein sequences
def fasta_sequence(fasta)
    first_seq = fasta.first.seq
    biosequence = Bio::Sequence.new(first_seq).auto   # guess the sequence type of this file
    if biosequence.instance_of? Bio::Sequence::NA
        return 'nucl'
    elsif biosequence.instance_of? Bio::Sequence::AA
        return 'prot'
    end
end

# Habría que hacer otra para comprobar si un archivo que sea tipo nucleótido es de CDSs

def only_cds?(nucleotide_fasta_file, start_codon, stop_codons)
    nucleotide_fasta_file.each do |entry|
        return false unless entry.seq.start_with?(start_codon)
        return false unless stop_codons.any? {|stop_codon| entry.seq.end_with?(stop_codon)}
    end
    return true
end
 
# Make database from proteome file depending of the type of sequences contained in the fasta file 

def make_db(fasta_file, filename, dbname)
    sequence_type = fasta_sequence(fasta_file)    # get the sequence type to define -dbtype
    result = system("makeblastdb -in #{filename} -dbtype '#{sequence_type}' -out ./Databases/#{dbname}")
    if $?.success?
        puts "BLAST database #{dbname} created successfully."
        puts
    else
        puts "Error creating BLAST database. Output:\n#{result}"
    end
end


# Ask user to input from comand line if wants to translate the CDS file for the proteome of some species

def ask_for_translation(file)
    puts "Do you want to translate #{file}? (yes/no):"
    response = gets.chomp.downcase
    case response
    when "no", "n"
        return false
    when "yes", "y"
        return true
    else 
        return "Invalid response. Please enter yes (y) or no (n)"
    end
end



# Translate CDS to protein: in case that user wants to perform reciprocal blastp

def search_sequence(identifier,filename)
    fasta = Bio::FlatFile.open(Bio::FastaFormat, filename)
    #puts "AAAA#{identifier}#{fasta}"
    fasta.each_entry do |entry|
        #puts "jajo #{entry.entry_id}.#{identifier}"
        #puts "UU #{entry.entry_id}"
        #MIRAR SI PONEEER EL IGUAAAAAL
        if entry.entry_id.match?(Regexp.escape(identifier)) #entry.definition.include?(identifier)
            #puts "hello #{identifier}"
          #puts "TIMEEEEEEEEEEEEEEEE #{i}"
          return entry
        end
    end
    puts "COULD NOT FIND YOUR SEQUENCE"
    return
end


def translate_cds_to_protein(cds_file, proteome_filename)
    proteome_file = File.open(proteome_filename, 'w')
  
    cds_file.each_entry do |entry|
      
      cds_seq = Bio::Sequence.auto(entry.seq)
      protein_seq = Bio::Sequence.auto(cds_seq.translate.chomp("*"))
      proteome_file.puts protein_seq.output_fasta(entry.definition) #entry.entry_id cuando busque
  
    end
      
end

# ------------------------ MAIN CODE --------------------------------------------

# Check input: introducing the file names of the proteomes of the species to discover putativo orthologues among 
if ARGV.length != 2
    abort "Incorrect number of files passed. Proteome files of the two species should be specified"
end

# Check order of arguments: to make the databases, we should know which file correspond to which species

if ARGV[0] == "TAIR10_cds.fa" && ARGV[1] == "proteome_pombe.fa"
    arabidopsis_cds = ARGV[0]
    pombe_proteome = ARGV[1]
elsif ARGV[1] == "TAIR10_cds.fa" && ARGV[0] == "proteome_pombe.fa"
    pombe_proteome = ARGV[0]
    arabidopsis_proteome = ARGV[1]
else
    abort "Incorrect files passed. Usage: ruby main.rb TAIR10_cds.fa proteome_pombe.fa"
end

arabidopsis_fasta = Bio::FlatFile.open(Bio::FastaFormat, arabidopsis_cds)

arabidopsis_proteome = "TAIR10_proteins.fa"
translate_cds_to_protein(arabidopsis_fasta, arabidopsis_proteome)

pombe_fasta = Bio::FlatFile.open(Bio::FastaFormat, pombe_proteome)
arabidopsis_fasta = Bio::FlatFile.open(Bio::FastaFormat, arabidopsis_proteome)

# 1st: Create both databases prior to performing reciprocal best BLAST
make_db(pombe_fasta, pombe_proteome, dbname = "POMBE")
make_db(arabidopsis_fasta, arabidopsis_proteome, dbname = "ARABIDOPSIS")

# 2nd: Prompt the user to specify from command line which type of search to do
# translate (y/n): if yes then a reciprocal blastp is performed, if not then tblastn + blastx


#if ask_for_translation(arabidopsis_proteome)
#    puts "ole"
#end

# 3rd: Create factories to perform blast depending on what kind of blast the user wants to perform

ara_factory = Bio::Blast.local('blastp', './Databases/ARABIDOPSIS')
pombe_factory = Bio::Blast.local('blastp', './Databases/POMBE')


# 4th: Perform blast and parse the output to do best hits reciprocal blast to find putative othologues

putative_othologues_candidates = Hash.new

EVALUE_THRESHOLD = 1e-6  # Threshold e-value as proposed in # meter ref que no me deja copiar y pegar en la maquina virtual


arabidopsis_fasta.each_entry do |entry|
    query = ">#{entry.entry_id}\n#{entry.seq}"
    report = pombe_factory.query(query)
  
    tblastn_best_hit = nil
  
    report.each_hit do |hit|
      if !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD
        tblastn_best_hit = hit if tblastn_best_hit.nil? || hit.evalue < tblastn_best_hit.evalue
        #puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
      end
    end
    #puts tblastn_best_hit
  
    next if tblastn_best_hit.nil?

    tblastn_best_hit_identifier = tblastn_best_hit.target_def.split("|")[0].delete(' ')
    #puts "aaaa #{tblastn_best_hit.query_def} #{tblastn_best_hit_identifier}"
    tblastn_best_hit_entry = search_sequence(tblastn_best_hit_identifier,pombe_proteome)
    query = ">#{tblastn_best_hit_entry.entry_id}\n#{tblastn_best_hit_entry.seq}"
    #puts
    #puts "#{tblastn_best_hit.target_def}#{tblastn_best_hit.target_seq}"
    #puts query
    #puts
    report = ara_factory.query(query)
  
    blastx_best_hit = nil
    report.each_hit do |hit|
      if !hit.evalue.nil? && hit.evalue <= EVALUE_THRESHOLD
  
        blastx_best_hit = hit if blastx_best_hit.nil? || hit.evalue < blastx_best_hit.evalue
        #puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
      end
    end
    #puts blastx_best_hit
  
    next if blastx_best_hit.nil?
    
    blastx_best_hit_identifier = blastx_best_hit.target_def.split("|")[0].delete(' ')

    if blastx_best_hit_identifier == entry.entry_id.delete(' ') # Si el blastx de la proteina->CDS da el ID del entry-CDS de TAIR
      puts "#{blastx_best_hit_identifier} is an orthologue candidate to #{tblastn_best_hit_entry.entry_id}"
      putative_othologues_candidates[blastx_best_hit_identifier] = tblastn_best_hit_identifier
    else
        puts "-#{blastx_best_hit_identifier}- is not equal to -#{entry.entry_id.delete(' ')}-"
    end
end


def write_candiates_report()

    File.open(output_report_file, 'w') do |file|

    file.puts "This code was created by Miguel La Iglesia Mirones and Lucía Muñoz Gil"
    file.puts "Bioinformatics Programming Challenges Course at MSc Computational Biology (UPM)"
    file.puts "December 2023"
    file.puts "-----------------------------------------------------------------------------------------------."
    file.puts
    file.puts "PUTATIVE ORTHOLOGUE CANDIDATES AMONG ARABIDOPSIS AND S.POMBE BY PERFORMING RECIPROCAL BEST BLAST" 
    file.puts
    file.puts "------------------------------------------------------------------------------------------------"
    file.puts
    file.puts
    file.puts "---------------------------------------------------------------------------------------"
    file.puts
    end
end 