from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import AlignIO

# Step 1: Perform BLAST Search
def perform_blast(sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    blast_record = NCBIXML.read(result_handle)
    # Extract top homologous sequences
    top_homologs = blast_record.alignments[:5]
    return top_homologs

# Step 2: Retrieve Homologous Sequences
def retrieve_sequences(top_homologs):
    homolog_sequences = []
    for alignment in top_homologs:
        for hsp in alignment.hsps:
            homolog_sequences.append(SeqRecord(seq=hsp.sbjct, id=alignment.title))
            break
    return homolog_sequences

# Step 3: Multiple Sequence Alignment
def perform_alignment(query_sequence, homolog_sequences):
    alignment_file = "alignment.fasta"
    SeqIO.write([query_sequence] + homolog_sequences, alignment_file, "fasta")
    clustalomega_cline = ClustalOmegaCommandline("clustalo", infile=alignment_file, outfile="alignment.clustal_num")
    stdout, stderr = clustalomega_cline()
    alignment = AlignIO.read("alignment.clustal_num", "clustal")
    return alignment

# Step 4: Calculate Conservation Score
def calculate_conservation_score(alignment):
    summary_align = AlignInfo.SummaryInfo(alignment)
    conservation_scores = summary_align.information_content_per_column()
    return conservation_scores

# Main function
def main():
    # Example sequence
    query_sequence = SeqIO.read("/Users/allalelhommad/PYT/SBPYT_project/dataset/protein_fasta/182L.fasta.txt", "fasta")

    # Step 1: Perform BLAST Search
    top_homologs = perform_blast(query_sequence.seq)

    print(top_homologs)
if __name__ == "__main__":
    main()
