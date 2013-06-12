annotation
==========
./fix_mRNA_lines.pl -f UNMC_rhesus_annotation_v2.12.gtf -o UNMC_rhesus_annotation_v2.13.gtf
./create_unique_tsids.pl UNMC_rhesus_annotation_v2.13.gtf 
gffread -y out -g MaSuRCA_Rhesus_Genome_v2.fasta UNMC_rhesus_annotation_v2.13.gtf.uniq 
./BLARGH/sort_gtf_file.pl UNMC_rhesus_annotation_v2.13.gtf.uniq
mv sorted.gtf UNMC_rhesus_annotation_v2.14.gtf
./pull_sequences_from_v2.pl -f MaSuRCA_Rhesus_Genome_v2.fasta -g UNMC_rhesus_annotation_v2.14.gtf
perl bowtie2.pl
perl shift_coordinates_with_bowtie2_results.pl -f UNMC_rhesus_annotation_v2.14.gtf -1 sam_R1 -2 sam_R2
perl ./BLARGH/sort_gtf_file.pl zz_new.gtf
perl ./BLARGH/fix_frames.pl sorted.gtf
perl ./fix_mRNA_lines.pl -f fixed_frames.gtf -o tmp.gtf
perl ./BLARGH/sort_gtf_file.pl tmp.gtf
gffread -y v3out -g MaSuRCA_Rhesus_Genome_v3_20130528.fasta sorted.gtf
gffread -y v3out_no_in-frame_stops -g MaSuRCA_Rhesus_Genome_v3_20130528.fasta sorted.gtf -V
