# NO FUNCIONA EN GENERAL LA BASE DE CDS DE TAIR
# MIRAR PARAMETROS
# MIRAR QUÉ HACER CUANDO TIENEN IGUAL E-VALUE



# REQUIRE PACKAGES
require  'bio'


# PUBLIC INSTANCE METHODS #
def make_blast_database(fasta_file,output_directory,database_name)
  makeblastdb_command = "makeblastdb -in #{fasta_file} -dbtype prot -out #{output_directory}/#{database_name}"
  result = `#{makeblastdb_command}`
  if $?.success?
    puts "BLAST database created successfully."
  else
    puts "Error creating BLAST database. Output:\n#{result}"
  end
end


# CONVERT CDS TO PROTEINS, I WON'T USE IT ...

def convert_cds_to_protein(cds_file)
  proteome_file = File.open('TAIR10_proteins.fa', 'w')

  cds_file.each_entry do |entry|
    
    cds_seq = Bio::Sequence.auto(entry.seq)
    protein_seq = Bio::Sequence.auto(cds_seq.translate.chomp("*"))
    proteome_file.puts protein_seq.output_fasta(entry.definition) #entry.entry_id cuando busque

  end
  proteome_file = Bio::FlatFile.auto('TAIR10_proteins.fa')
  return proteome_file

end

#cds_file = Bio::FlatFile.auto('TAIR10_cds.fa')
#tair_file = convert_cds_to_protein(cds_file)


# MAIN CODE #

tair_filename = 'TAIR10_cds.fa'
pep_filename = 'pep.fa'

tair_file = Bio::FlatFile.auto(tair_filename)
pep_file = Bio::FlatFile.auto(pep_filename)


make_blast_database(tair_filename,'databases','TAIR10')
make_blast_database(pep_filename,'databases','pep')


tair_database = Bio::Blast.local('tblastn', './databases/TAIR10')
pep_database = Bio::Blast.local('blastx', './databases/pep')

evalue_threshold = 1e-10  # Adjust this value as needed

tair_file.each_entry do |entry|
  query = ">#{entry.entry_id}\n#{entry.seq}"
  report = pep_database.query(query)

  tblastn_best_hit = nil

  report.each_hit do |hit|
    if !hit.evalue.nil? && hit.evalue <= evalue_threshold
      tblastn_best_hit = hit if tblastn_best_hit.nil? || hit.evalue < tblastn_best_hit.evalue
      #puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
    end
  end
  #puts tblastn_best_hit

  next if tblastn_best_hit.nil?

  query = ">#{tblastn_best_hit.hit_id}\n#{tblastn_best_hit.target_seq}"
  puts query
  report = tair_database.query(query)

  blastx_best_hit = nil

  report.each_hit do |hit|
    if !hit.evalue.nil? && hit.evalue <= evalue_threshold

      blastx_best_hit = hit if blastx_best_hit.nil? || hit.evalue < blastx_best_hit.evalue
      #puts "#{hit.hit_id} : evalue #{hit.evalue}\t#{hit.target_id} "
    end
  end
  #puts blastx_best_hit

  next if blastx_best_hit.nil?

  if blastx_best_hit.hit_id == entry.entry_id # Si el blastx de la proteina->CDS da el ID del entry-CDS de TAIR
    puts "#{entry.entry_id} is orthologue to #{tblastn_best_hit.hit_id}"
  end

end