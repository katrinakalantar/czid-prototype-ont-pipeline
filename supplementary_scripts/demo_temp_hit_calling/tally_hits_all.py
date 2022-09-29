#
# example of how to run this (specifically in the demo folder):
# python3 tally_hits_all.py 20 20220823_145511_ontpipeline.reads-to-contigs.out
#

import sys
import pandas as pd

file_root = sys.argv[1]
m8_df = pd.read_csv(file_root + ".gsnap.deduped.m8", sep="\t", names = ["read_id", "accession", "pid", "aln_len", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"])
hitsummary_df = pd.read_csv(file_root + ".gsnap.hitsummary.tab", sep="\t", names = ["read_id", "x", "final_taxid", "accession", "species_taxid", "genus_taxid", "family_taxid"])

big_df = pd.concat([m8_df, hitsummary_df], axis=1)
sub_df = big_df[["read_id","aln_len","final_taxid", "genus_taxid"]]
sub_df.columns = ["read_id", "read_id2", "aln_len","final_taxid", "genus_taxid"]

# maybe this should be an ASSERT statement
if not (sum(sub_df["read_id"] == sub_df["read_id2"]) == len(sub_df.index)):
	print("ERROR: sequence IDs don't match")

sp_data_collection = sub_df[["final_taxid","aln_len"]]
sp_counts = sp_data_collection.groupby(["final_taxid"]).sum()

# this should contain the sorted species counts based on reads
print(sp_counts.sort_values(by="aln_len", ascending = False).head(10))

gen_data_collection = sub_df[["genus_taxid","aln_len"]]
gen_counts = gen_data_collection.groupby(["genus_taxid"]).sum()
# this should contain the sorted genus counts based on reads
print(gen_counts.sort_values(by="aln_len", ascending = False).head(10))


# note: this file comes from modification of the RunReadsToContigs .sam file
# less 20220722_022439_ontpipeline/call-RunReadsToContigs/work/sample.reads_to_contigs.sam | grep -v "^@" | cut -f1,3,10 > 20220722_022439_ontpipeline.reads-to-contigs.out
# so, we would need to take in the the reads_to_contigs.sam and remove the header and then take just those columns (1,3,10) to create the input for this section
r2c_df = pd.read_csv(sys.argv[2], 
                     sep='\t', names = ['seqid', 'contig', 'sequence'])
r2c_df['seq_len'] = [len(i) for i in r2c_df['sequence']]
r2c_df = r2c_df[r2c_df.contig != '*']
r2c_df.head()
contig_tally_df = r2c_df[['contig', 'seq_len']]
result = contig_tally_df.groupby(['contig']).sum()

sub_df.join(result)
df = pd.merge(sub_df, result, how = 'outer', left_on="read_id", right_on="contig")
sp_df = df[['final_taxid', 'seq_len']]
gen_df = df[['genus_taxid', 'seq_len']]

# these should contain the sorted species / genus counts by length (bp)
sp_result_final = sp_df.groupby(["final_taxid"]).sum()
gen_result_final = gen_df.groupby(["genus_taxid"]).sum()

print(sp_result_final.sort_values(by="seq_len", ascending = False).head())
print(gen_result_final.sort_values(by="seq_len", ascending = False).head())

sp_counts.to_csv("output.sp_counts.csv")
gen_counts.to_csv("output.gen_counts.csv")
sp_result_final.to_csv("output.sp_result_final.csv")
gen_result_final.to_csv("output.gen_result_final.csv")