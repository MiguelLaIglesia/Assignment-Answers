class EmblEntry

    attr_accessor :entry_sequence, :chromosomal_coordinates, :source, :accession, :seq_id, :gene_locus

    def initialize(entry, start_entry, end_entry)
        @accession = entry.accession
        @seq_id = entry.entry_id
        @entry_sequence = entry.to_biosequence
        @chromosomal_coordinates = [start_entry, end_entry]
        @source = ""
        @gene_locus =""
    end

    def annotate_source_gene
        @entry_sequence.features.each do |feature|
          featuretype = feature.feature
          position = feature.position
          qual = feature.assoc
          
          # Annotate source
          if featuretype == 'source'
            @source = qual['db_xref']
          end
          
          # Annotate gene locus for the sequence
          if featuretype == 'gene' && position.match(/^(complement\()?(\d+\.\.\d+)(\))?$/) && @gene_locus.empty?
            @gene_locus = qual['gene']
          end

        end
        return
    end


    # Add found repeats as feature of the Embl biosequence object
    def add_cttctt_repeat(coordinates, start_motif, end_motif, strand, source)
        return if @entry_sequence.features.any? { |f| f.feature == 'cttctt_repeat' && f.position == coordinates }   # check if this motif has already been annotated
    
        f1 = Bio::Feature.new('cttctt_repeat', coordinates) # New featyre of the Embl biosequence object --> repetitive motifs "CTTCTT"
        f1.append(Bio::Feature::Qualifier.new('start', start_motif))
        f1.append(Bio::Feature::Qualifier.new('end', end_motif))
        f1.append(Bio::Feature::Qualifier.new('strand', strand))
        f1.append(Bio::Feature::Qualifier.new('SO_Name', 'repeat_region'))
        @entry_sequence.features << f1
    end

    # convert motif sequnce-coordenates into chromosomal-coordenates for second part of the task
    def convert_to_chromosomal_coordinates(start_motif, end_motif)
        start_chr = @chromosomal_coordinates[0]
        length_motif = end_motif - start_motif + 1
        start_motif = start_motif + start_chr - 1
        end_motif = start_motif + length_motif - 1
        return [start_motif, end_motif]
    end

    def process_cttctt_repeats(coordinates_to_use, regex_positive=REGEX_POSITIVE, regex_complementary = REGEX_COMPLEMENT)

        @entry_sequence.features.each do |feature|
          next unless feature.feature == 'exon' && feature.position.match?(/^(complement\()?(\d+\.\.\d+)(\))?$/)    # check for exons
    
          start_exon, end_exon = get_coordinates(feature.position)
          exon_sequence = @entry_sequence.subseq(start_exon, end_exon)  # get exon subsequence from exon coordinates
    
          next if exon_sequence.nil? || exon_sequence.empty?    # control: is empty?
    
          if feature.position.include?('complement')    # get directionality of exon (positive or complement strand)
            regex = regex_complementary
            strand = '-'
            coordinates_format = "complement(%s)"
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
            start_motif = start_exon + start_match  # +1 por ccs del match, -1 por la suma de exon y match
            end_motif = start_exon + end_match - 1  # -1 por la suma de exon y match
    
             # If chromosomal coordinates, it recalculates those motif coordinates in the chromosomal coordinates
            if coordinates_to_use == 'chromosomal coordinates'
                start_motif, end_motif = convert_to_chromosomal_coordinates(start_motif, end_motif)
            end

            coordinates = coordinates_format % "#{start_motif}..#{end_motif}"
            add_cttctt_repeat(coordinates, start_motif, end_motif, strand, @source)
            
          end
        end
      end

end