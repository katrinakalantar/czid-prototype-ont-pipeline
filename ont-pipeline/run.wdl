version 1.0

workflow phylotree{
        input{
                File input_fastq
                String library_type = "RNA"
                String guppy_basecaller_version = "5.1"
                String guppy_basecaller_setting = "hac" # fast, hac, super

                # only necessary if using medaka for polishing contigs
                #String platform = "minION" # minION, gridION, promethION
                #String pore = "r9.4.1"

                String subsample_depth = 4000000 # should be 4x the number of reads desired

                File minimap_host_db
                File minimap_human_db

                File NT_minimap2
                File NT_centrifuge
                String alignment_test_mode = "all_mm" # all_mm, split_mm_cent 

                File NR_diamond
                
                String polishing_iterations = 1

                String docker_image_id
        }

        call RunValidateInput{
                input:
                        input_fastq = input_fastq,
                        docker_image_id = docker_image_id
        }

        call RunQualityFilter{
                input:
                        input_fastq = input_fastq,
                        docker_image_id = docker_image_id
        }

        call RunHostFilter{
                input:
                        input_fastq = RunQualityFilter.fastp_output,
                        library_type = library_type,
                        minimap_host_db = minimap_host_db,
                        docker_image_id = docker_image_id
        }

        call RunHumanFilter{
                input:
                        input_fastq = RunHostFilter.host_filter_fastq,
                        library_type = library_type,
                        minimap_human_db = minimap_human_db,
                        docker_image_id = docker_image_id
        }

        call RunSubsampling{
                input:
                        input_fastq = RunHumanFilter.human_filter_fastq,
                        subsample_depth = subsample_depth,
                        docker_image_id = docker_image_id
        }

        call RunAssembly{
                input:
                        input_fastq = RunSubsampling.subsampled_fastq,
                        guppy_basecaller_setting = guppy_basecaller_setting,
                        polishing_iterations = polishing_iterations,
                        docker_image_id = docker_image_id
        }

        call RunReadsToContigs{
                input:
                        input_fastq = RunSubsampling.subsampled_fastq,
                        assembled_reads = RunAssembly.assembled_fasta,
                        docker_image_id = docker_image_id
        }

        call RunNTAlignment{
                input:
                        non_contig_reads_fa = RunReadsToContigs.non_contigs_fasta,
                        assembled_reads_fa = RunAssembly.assembled_fasta,
                        NT_minimap2 = NT_minimap2,
                        NT_centrifuge = NT_centrifuge,
                        alignment_test_mode = alignment_test_mode,
                        docker_image_id = docker_image_id
        }

        call RunNRAlignment{
                input:
                        assembled_reads_fa = RunAssembly.assembled_fasta,
                        NR_diamond = NR_diamond,
                        docker_image_id = docker_image_id
        }

        output{
               File validate_input_dummy_output = RunValidateInput.dummy_output
               File validate_input_fastq = RunValidateInput.validated_output
               File dummy_output = RunQualityFilter.dummy_output
               File fastp_output = RunQualityFilter.fastp_output
               File host_dummy_output = RunHostFilter.dummy_output
               File hostfilter_fastq = RunHostFilter.host_filter_fastq
               File hostfilter_sam = RunHostFilter.host_filter_sam
               File human_dummy_output = RunHumanFilter.dummy_output
               File humanfilter_fastq = RunHumanFilter.human_filter_fastq
               File humanfilter_sam = RunHumanFilter.human_filter_sam
               File subsample_dummy_output = RunSubsampling.dummy_output
               File subsampled_fastq = RunSubsampling.subsampled_fastq
               File assembled_fasta = RunAssembly.assembled_fasta
               File assembly_dir = RunAssembly.temp_assembly_dir
               File reads_to_contigs_dummy_output = RunReadsToContigs.dummy_output
               File reads_to_contigs_non_contigs_fastq = RunReadsToContigs.non_contigs_fastq
               File reads_to_contigs_non_contigs_fasta = RunReadsToContigs.non_contigs_fasta
               File run_nt_alignment_dummy_output = RunNTAlignment.dummy_output
               File run_nt_alignment_all_seqs = RunNTAlignment.all_sequences_to_align
               File run_nt_minimap_sam = RunNTAlignment.nt_minimap_sam
               File run_nt_centrifuge_output = RunNTAlignment.nt_centrifuge_output
               File run_nr_alignment_dummy_output = RunNRAlignment.dummy_output
               File run_nr_alignment_diamond_output = RunNRAlignment.diamond_output
        }
}

task RunValidateInput{
        input{
               File input_fastq
               String docker_image_id
        }
        command <<<
               echo "inside RunValidateInput step" >> output.txt
               # TODO: add some validation, in the meantime just copy to output
               cp "~{input_fastq}" sample.validated
        >>>

        output{
               File dummy_output = "output.txt"
               File validated_output = "sample.validated"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunQualityFilter{
        input{
               File input_fastq
               String docker_image_id
        }
        command <<<
               echo "inside RunQualityFilter step" >> output.txt
               fastp -i "~{input_fastq}" --qualified_quality_phred 9 --length_required 100 --low_complexity_filter --complexity_threshold 30 --dont_eval_duplication -o sample.fastp
        >>>

        output{
               File dummy_output = "output.txt"
               File fastp_output = "sample.fastp"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunHostFilter{
        input{
               File input_fastq
               String library_type
               File minimap_host_db
               String docker_image_id
        }
        command <<<
               echo "inside RunHostFilter step" >> output.txt
               # run minimap2 against host genome
               minimap2 -ax splice "~{minimap_host_db}" "~{input_fastq}" -o sample.hostfiltered.sam
               # extract the unmapped reads for downstream processing
               samtools fastq -n -f 4 sample.hostfiltered.sam > sample.hostfiltered.fastq
        >>>

        output{
               File dummy_output = "output.txt"
               File host_filter_sam = "sample.hostfiltered.sam"
               File host_filter_fastq = "sample.hostfiltered.fastq"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunHumanFilter{
        input{
               File input_fastq
               String library_type
               File minimap_human_db
               String docker_image_id
        }
        command <<<
               echo "inside RunHumanFilter step" >> output.txt
               # run minimap2 against human genome
               minimap2 -ax splice "~{minimap_human_db}" "~{input_fastq}" -o sample.humanfiltered.sam
               # extract the unmapped reads for downstream processing
               samtools fastq -n -f 4 sample.humanfiltered.sam > sample.humanfiltered.fastq
        >>>

        output{
               File dummy_output = "output.txt"
               File human_filter_sam = "sample.humanfiltered.sam"
               File human_filter_fastq = "sample.humanfiltered.fastq"
        }
        runtime{
               docker: docker_image_id
        }
}


task RunSubsampling{
        input{
               File input_fastq
               String subsample_depth
               String docker_image_id
        }
        command <<<
               echo "inside RunSubsampling step" >> output.txt
               head -"~{subsample_depth}" "~{input_fastq}" > sample.subsampled.fastq
        >>>

        output{
               File dummy_output = "output.txt"
               File subsampled_fastq = "sample.subsampled.fastq"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunAssembly{
        input{
               File input_fastq
               String guppy_basecaller_setting
               String polishing_iterations
               String docker_image_id
        }
        command <<<
               echo "inside RunAssembly step" >> output.txt
               flye_setting="--nano-raw"
               if ["~{guppy_basecaller_setting}" -eq "super"]
               then
                 echo "DEBUG: inside loop for flye_setting is set to --nano-hq" >> output.txt
                 flye_setting="--nano-hq"
               fi

               echo "DEBUG: flye_setting = $flye_setting" >> output.txt

               # run flye to assembly contigs
               flye --meta $flye_setting "~{input_fastq}" --out-dir temp_flye_out --threads 8 --iterations "~{polishing_iterations}"

               # ERROR HANDLING - assembly somethings fails (due to low coverage) and is then missing...
               #                  ... the temp_flye_out/assembly.fasta file
               if [ -f temp_flye_out/assembly.fasta ]
               then
                 echo "DEBUG: assembly.fasta file is found" >> output.txt
                 cat temp_flye_out/assembly.fasta > sample.assembled_reads.fasta
               else
                 echo "DEBUG: assembly.fasta file is not found, copying input to output" >> output.txt
                 #just copy original .fastq to .fasta
                 seqtk seq -a "~{input_fastq}" > sample.assembled_reads.fasta 
               fi

               zip -r temp_flye_out.zip temp_flye_out

        >>>

        output{
               File dummy_output = "output.txt"
               File assembled_fasta = "sample.assembled_reads.fasta"
               File temp_assembly_dir = "temp_flye_out.zip"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunReadsToContigs{
        input{
               File input_fastq
               File assembled_reads
               String docker_image_id
        }
        command <<<
               echo "inside RunReadsToContigs step" >> output.txt

               # use minimap2 to align reads back to contigs
               minimap2 -ax map-ont "~{assembled_reads}" "~{input_fastq}" -o sample.reads_to_contigs.sam

               # extract reads that didn't map to contigs as non-contig reads 
               samtools fastq -n -f 4 sample.reads_to_contigs.sam > sample.non_contigs.fastq

               # convert non-contigs.fastq file to .fasta file
               seqtk seq -a sample.non_contigs.fastq > sample.non_contigs.fasta

        >>>

        output{
               File dummy_output = "output.txt"
               File non_contigs_fastq = "sample.non_contigs.fastq"
               File non_contigs_fasta = "sample.non_contigs.fasta"
        }
        runtime{
               docker: docker_image_id
        }
}


task RunNTAlignment{
        input{
               File non_contig_reads_fa
               File assembled_reads_fa
               File NT_minimap2
               File NT_centrifuge
               String alignment_test_mode
               String docker_image_id
        }
        command <<<
               echo "inside RunNTAlignment step" >> output.txt

               echo "seqs in non_contigs.fasta" >> output.txt
               grep ">" "~{non_contig_reads_fa}" | wc -l >> output.txt

               # combine contigs and non-contig reads into a single file
               cat "~{assembled_reads_fa}" "~{non_contig_reads_fa}" > sample.all_sequences_to_align_full.fasta

               echo "seqs in all_sequences_to_align_full.fasta" >> output.txt
               grep ">" sample.all_sequences_to_align_full.fasta | wc -l >> output.txt

               # remove duplicates from full .fasta (important when assembly failed)
               seqkit rmdup -s < sample.all_sequences_to_align_full.fasta > sample.all_sequences_to_align.fasta

               echo "seqs in all_sequences_to_align.fasta" >> output.txt
               grep ">" sample.all_sequences_to_align.fasta | wc -l >> output.txt

               alignment_mode="~{alignment_test_mode}"
               echo "alignment mode = $alignment_mode" >> output.txt

               # IF running option #1, map contigs and non-contig reads to NT with minimap
               if [[ $alignment_mode == "all_mm" ]]
               then
                 echo "DEBUG: alignment mode = $alignment_mode; inside full minimap alignment" >> output.txt
                 minimap2 -ax asm20 -o sample.nt_minimap2_output.sam "~{NT_minimap2}" sample.all_sequences_to_align.fasta
                 # create dummy centrifuge output 
                 touch sample.nt_centrifuge_output.txt
               # IF running option #2...
               # map assembled_reads_fa to NT with minimap and map non_contig_reads_fa to NT with centrifuge
               elif [[ $alignment_mode == "split_mm_cent" ]]
               then
                 echo "DEBUG: alignment mode = $alignment_mode; inside split alignment" >> output.txt
                 # Run minimap2 on just the contigs
                 minimap2 -ax asm20 -o sample.nt_minimap2_output.sam "~{NT_minimap2}" "~{assembled_reads_fa}"                 
                 # Run centrifuge on non-contig reads
                 unzip "~{NT_centrifuge}"
                 centrifuge_dir=`basename -s .zip reference/centrifuge-ref.zip`
                 centrifuge -f --min-hitlen 50 -x "reference/$centrifuge_dir/p_compressed+h+v" -U "~{non_contig_reads_fa}" -S sample.nt_centrifuge_output.txt
               else
                 echo "ERROR: UNKNOWN VALUE FOR alignment_test_mode" >> output.txt
               fi
        >>>

        output{
               File dummy_output = "output.txt"
               File all_sequences_to_align = "sample.all_sequences_to_align.fasta"
               File nt_minimap_sam = "sample.nt_minimap2_output.sam"
               File nt_centrifuge_output = "sample.nt_centrifuge_output.txt"
        }
        runtime{
               docker: docker_image_id
        }
}

task RunNRAlignment{
        input{
               File assembled_reads_fa
               File NR_diamond
               String docker_image_id
        }
        command <<<
               echo "inside RunNRAlignment step" >> output.txt
               # run DIAMOND
               diamond blastx --long-reads -d "~{NR_diamond}" -q "~{assembled_reads_fa}" -o sample.nr_diamond_output.txt
        >>>

        output{
               File dummy_output = "output.txt"
               File diamond_output = "sample.nr_diamond_output.txt"
        }
        runtime{
               docker: docker_image_id
        }
}
