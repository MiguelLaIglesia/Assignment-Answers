# == EmblProcessor
#
# This is a representation of an EMBL file processor
# which scans entry features for creating GFF file
#
# == Summary
# 
# This can be used to represent embl files which will be processed for creating GFF with embl entries features
#

class EmblProcessor
    
    # Get/Set total entries in embl file
    # @!attribute [rw]
    # @return [String] all entries
    attr_accessor :entries

    
    # Create a new instance of EmblProcessor
    # @param filename [String] the name of the embl file
    # @return [EmblProcessor] an instance of EmblProcessor
    def initialize(filename)
        @entries = scan_entries(filename)
    end

    # Scan embl file entries to create  EmblEntry instance and save information about it
    # @param filename [String] the name of the embl file
    # @return [entries] all entries in embl file
    def scan_entries(filename)
        embl_file = Bio::FlatFile.auto(filename)
        entries = []

        embl_file.each_entry do |entry|
            next unless entry.accession # Lack of accesion is suspicious
            start_entry, end_entry = get_coordinates(entry.definition)  #  Extract sequence chromosome coordinates  
            embl_entry = EmblEntry.new(entry, start_entry, end_entry)
            embl_entry.annotate_source_gene # Annotates source and gene for entry
            entries << embl_entry   # All entries in embl file
        end
        
        return entries
    end  

    # Creates gff for embl file with repetitive regions as features
    # @param coordinates_to_use [String] the coordinates of reference to consider
    # @return [gff] the GFF file
    def load_to_gff(coordinates_to_use)
        
        report_count = 0 # Number of entries without 'cttctt'
        entries_without_regions = []

        gff = Bio::GFF::GFF3.new
        feature_count = 1
    
        @entries.each do |embl_entry|

          # Creates new feature with repetitive regions
          embl_entry.process_cttctt_repeats(coordinates_to_use) 
          
          # Test if entry has repetitive regions annotated
          feature_exists = embl_entry.entry_sequence.features.any? { |feature| feature.feature == 'cttctt_repeat' }
          
          if !feature_exists
            entries_without_regions << embl_entry
            report_count += 1
            next
          end

          embl_entry.entry_sequence.features.each do |feature|

            next unless feature.feature == 'cttctt_repeat'
    
            qual = feature.assoc
            attributes = [{ 'ID' => "repeat_region_#{feature_count}", 'Name' => "cttctt_repeat_#{embl_entry.gene_locus}" }]
    
            seqid = embl_entry.seq_id
            
            # New record for GFF file
            gff.records << Bio::GFF::GFF3::Record.new(
              seqid,            # seqID
              'programatically',# source
              qual['SO_Name'],  # feature type
              qual['start'],    # start
              qual['end'],      # end
              nil,              # score
              qual['strand'],   # strand
              nil,              # phase
              attributes[0]     # attributes
            )
            feature_count += 1
          end
        end

        # Writes a report with genes without repetitive regions
        report = File.open("gff_report.txt", "w")
        report.puts "REPORT OF GENES WITHOUT \'CTTCTT\' REPETITIVE REGIONS AFTER CREATING GFF FILE"
        report.puts "---------------------------------------------------------------------------"
        entries_without_regions.each do |entry|
          report.puts "#{entry.gene_locus}"
        end
        report.puts
        report.puts "TOTAL GENES: #{report_count}"
        report.close

        return gff
      end


    
end