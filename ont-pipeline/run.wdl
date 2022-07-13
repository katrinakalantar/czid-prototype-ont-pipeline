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
                
                String polishing_iterations = 1

                String docker_image_id
        }

        call RunHello{
                input:
                        input_fastq = input_fastq,
                        docker_image_id = docker_image_id
        }

        call RunHostFilter{
                input:
                        input_fastq = RunHello.fastp_output,
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

        output{
               File dummy_output = RunHello.dummy_output
               File fastp_output = RunHello.fastp_output
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
        }
}

task RunHello{
        input{
               File input_fastq
               String docker_image_id
        }
        command <<<
               echo "hello world" >> output.txt
               centrifuge --help >> output.txt
               head "~{input_fastq}" >> output.txt
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
               echo "hello world from inside host filter" >> output.txt
               minimap2 -ax splice "~{minimap_host_db}" "~{input_fastq}" -o sample.hostfiltered.sam
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
               echo "hello world from inside host filter" >> output.txt
               minimap2 -ax splice "~{minimap_human_db}" "~{input_fastq}" -o sample.humanfiltered.sam
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
               echo "hello world from inside host filter" >> output.txt
               head "~{subsample_depth}" "~{input_fastq}" > sample.subsampled.fastq
               echo "lines in input:" >> output.txt
               wc -l "~{input_fastq}" >> output.txt
               echo "lines in output:" >> output.txt
               wc -l sample.subsampled.fastq >> output.txt
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
               echo "hello world from inside assembly step" >> output.txt
               flye_setting="--nano-raw"
               if ["~{guppy_basecaller_setting}" -eq "super"]
               then
                 echo "inside loop" >> output.txt
                 flye_setting="--nano-hq"
               fi

               echo "flye_setting = " >> output.txt
               echo $flye_setting >> output.txt
               echo "$flye_setting" >> output.txt

               flye --help >> output.txt

               flye --meta $flye_setting "~{input_fastq}" --out-dir temp_flye_out --threads 8 --iterations "~{polishing_iterations}"

               # ERROR HANDLING - ASSEMBLY SOMETIMES FAILS TO COMPLETE AND IS THEN MISSING THE temp_flye_out/assembly.fasta file
               if [ -f temp_flye_out/assembly.fasta ]
               then
                 echo "File is found" >> output.txt
                 cat temp_flye_out/assembly.fasta > sample.assembled_reads.fasta
               else
                 echo "File is not found" >> output.txt
                 seqtk seq -a "~{input_fastq}" > sample.assembled_reads.fasta #just copy original .fastq to .fasta
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
               echo "hello world from inside ReadsToContigs" >> output.txt

               echo "seqs in input reads .fastq" >> output.txt
               grep "^@" "~{input_fastq}" | wc -l >> output.txt


               minimap2 -ax map-ont "~{assembled_reads}" "~{input_fastq}" -o sample.reads_to_contigs.sam
               samtools fastq -n -f 4 sample.reads_to_contigs.sam > sample.non_contigs.fastq

               echo "seqs in assembed reads .fasta" >> output.txt
               grep ">" "~{assembled_reads}" | wc -l >> output.txt
               
               echo "seqs in non_contigs.fastq" >> output.txt
               grep "^@" sample.non_contigs.fastq | wc -l >> output.txt

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
               echo "hello world from inside RunNTAlignment" >> output.txt

               echo "seqs in non_contigs.fasta" >> output.txt
               grep ">" "~{non_contig_reads_fa}" | wc -l >> output.txt

               # COMBINE CONTIGS and NON-CONTIG READS INTO file
               cat "~{assembled_reads_fa}" "~{non_contig_reads_fa}" > sample.all_sequences_to_align_full.fasta

               echo "seqs in all_sequences_to_align_full.fasta" >> output.txt
               grep ">" sample.all_sequences_to_align_full.fasta | wc -l >> output.txt

               # remove duplicates
               seqkit rmdup -s < sample.all_sequences_to_align_full.fasta > sample.all_sequences_to_align.fasta

               echo "seqs in all_sequences_to_align.fasta" >> output.txt
               grep ">" sample.all_sequences_to_align.fasta | wc -l >> output.txt

               echo "alignment test mode is: " >> output.txt
               echo "~{alignment_test_mode}" >> output.txt
               alignment_mode="~{alignment_test_mode}"
               echo $alignment_mode >> output.txt
               echo $(( $alignment_mode == "all_mm" )) >> output.txt


               if $(( $alignment_mode == "all_mm" ))
               then
                 echo "HELLOOOOO - test1" >> output.txt
               fi

               if [[ $alignment_mode == "all_mm" ]]
               then
                 echo "HELLOOOOO - test2" >> output.txt
               fi

               if [[ $alignment_mode == "all_mm" ]]
               then
                 echo "HELLOOOOO" >> output.txt
               elif [[ $alignment_mode == "split_mm_cent" ]]
               then
                 echo "HELLLOOOOO2" >> output.txt
               fi



               # IF RUNNING OPTION 1, map sample.all_sequences_to_align.fasta to NT with minimap
               if [[ $alignment_mode == "all_mm" ]]
               then
                 echo "inside full minimap alignment" >> output.txt
                 minimap2 -ax asm20 -o sample.nt_minimap2_output.sam "~{NT_minimap2}" sample.all_sequences_to_align.fasta
                 touch sample.nt_centrifuge_output.txt
               # IF RUNNING OPTION 2, map assembled_reads_fa to NT with minimap and map non_contig_reads_fa to NT with centrifuge
               elif [[ $alignment_mode == "split_mm_cent" ]]
               then
                 echo "inside split alignment of minimap and centrifuge" >> output.txt
                 # Run minimap2 on just the contigs
                 minimap2 -ax asm20 -o sample.nt_minimap2_output.sam "~{NT_minimap2}" "~{assembled_reads_fa}"
                 # Run centrifuge on non-contig reads

                 unzip "~{NT_centrifuge}"
                 centrifuge_dir=`basename -s .zip reference/centrifuge-ref.zip`
                 echo $centrifuge_dir >> output.txt
                 echo "reference/$centrifuge_dir/p_compressed+h+v" >> output.txt
                 ls $centrifuge_dir >> output.txt

                 centrifuge -q --min-hitlen 50 -x "reference/$centrifuge_dir/p_compressed+h+v" -U "~{non_contig_reads_fa}" -S sample.nt_centrifuge_output.txt
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


