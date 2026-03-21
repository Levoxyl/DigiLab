import os 
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def find_virus_parts(input_fasta, target_root):
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    output_dir = os.path.join(target_root, f"{base_name}_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        all_proteins = str(record.seq.translate(to_stop=False)).split('*')
        real_proteins = [p for p in all_proteins if len(p) > 0]
        
        for i, p_seq in enumerate(real_proteins):
            analysed_p = ProteinAnalysis(p_seq)
            
            weight = analysed_p.molecular_weight()
            gravy = analysed_p.gravy()
            
            location = "External (Water-Loving)" if gravy < 0 else "Internal/Membrane (Water-Hating)"    
            
            file_path = os.path.join(output_dir, f"part_{i}.txt")
            with open(file_path, 'w') as f:
                f.write(f"Source: {record.description}\n")
                f.write(f"Fragment index: {i}\n")
                f.write(f"Molecular Weight: {weight:.2f} Daltons\n")
                f.write(f"GRAVY Score: {gravy:.2f}\n")
                f.write(f"Location Prediction: {location}\n")
                f.write(f"-" * 20 + "\n")
                f.write(p_seq)
                
    print(f"Done! Created {len(real_proteins)} files in {output_dir}/")
