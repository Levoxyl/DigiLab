import os
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def process_lab_directory(input_fasta):
    base_name = os.path.splitext(input_fasta)[0]
    organism_folder = f"{base_name}_proteins"
    os.makedirs(organism_folder, exist_ok=True)

    print(f"Starting analysis for: {input_fasta}")

    for record in SeqIO.parse(input_fasta, "fasta"):
        full_protein = record.seq.translate()
        all_proteins = str(full_protein).split('*')
        real_proteins = [p for p in all_proteins if len(p) > 20]

        for i , p_seq in enumerate(real_proteins):
            try:
                analysed_p = ProteinAnalysis(p_seq)
                weight = analysed_p.molecular_weight()
                file_name = f"{base_name}_fragment_{i}.txt"
                file_path = os.path.join(organism_folder, file_name)

                with open(file_path, 'w') as f:
                    f.write(f"Source: {record.description}\n")
                    f.write(f"Fragment index: {i}\n")
                    f.write(f"Molecular Weight: {weight:.2f} Daltons\n")
                    f.write(f"-" * 20 + "\n")
                    f.write(p_seq)
            except Exception as e:
                print(f"Skipped fragment {i} due to calculation error.")

            print(f"Done! Created {len(real_proteins)} files in {organism_folder}/")

process_lab_directory("SARS_CoV_2.fasta")

# Partial codon warning
# If DNA sequence length is not a multiple of 3 -> the translate function  
# Will issue a warning about a partial codon at the end of the sequence.
# The translation process reads the DNA in groups of 3 nucleotides (codons) to produce amino acids
# If there are leftover nucleotides that do not form a complete codon
# It can lead to an incomplete translation. 
#
# To avoid this warning
# You can ensure that your DNA sequence length is a multiple of 3 
# Or 
# Handle the partial codon appropriately in your code.