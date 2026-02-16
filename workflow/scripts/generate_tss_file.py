import sys
import pyranges1 as pr

"""
Note to self about 1-based GTF vs 0-based BED:

Some genomic coordinate formats (GFF, GTF) use a 1-based, start-and-end-included format. Pyranges takes care of converting between these conventions when loading and writing files in these formats.

Source:
https://pyranges1.readthedocs.io/en/latest/tutorial.html

"""

# Get command line arguments
args = sys.argv[1:]
if len(args) != 2:
    print("Usage: python generate_tss_file.py <annotations.gtf> <output.bed>")
    sys.exit(1)

# Load the GTF
gr = pr.read_gtf(args[0])

# Filter for transcripts
transcripts = gr[gr.Feature == "transcript"]

# Remove mitochondrial transcripts
transcripts = transcripts[transcripts.Chromosome != "MT"]

# Calculate TSS
# .five_end() handles the strand-specific logic
tss = transcripts.five_end()

# Combine Gene ID and Transcript ID
tss["Name"] = tss.gene_id + "|" + tss.transcript_id

# Add a dummy score
tss["Score"] = 0

# Select only the columns needed for BED
tss = tss[["Chromosome", "Start", "End", "Name", "Score", "Strand"]]

# Export to BED
# PyRanges automatically maps 'Name' to the 4th column of the BED file
output_bed = args[1]
tss.to_bed(output_bed)
