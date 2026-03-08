import tkinter as tk
from tkinter import filedialog, messagebox
import os
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class BioWorkbench:
    def __init__(self, root):
        self.root = root
        self.root.title("GENOMIC_ANALYZER_V1.1")
        self.root.geometry("600x550")
        self.root.configure(bg="#d9d9d9") # Classic Windows 95 Grey

        self.fasta_path = tk.StringVar()
        self.output_dir = tk.StringVar()

        # --- HEADER ---
        tk.Label(root, text="LAB_TERMINAL: DNA_DIGEST_SYSTEM", bg="#000080", fg="white", 
                 font=("Courier", 12, "bold"), relief="raised", bd=3).pack(fill="x", pady=5)

        # --- INPUTS ---
        input_frame = tk.Frame(root, bg="#d9d9d9", bd=2, relief="groove")
        input_frame.pack(fill="x", padx=10, pady=10)

        tk.Button(input_frame, text="FILE...", command=self.select_file, width=10).grid(row=0, column=0, padx=5, pady=5)
        tk.Entry(input_frame, textvariable=self.fasta_path, width=50, bg="white").grid(row=0, column=1, padx=5)

        tk.Button(input_frame, text="DIR...", command=self.select_folder, width=10).grid(row=1, column=0, padx=5, pady=5)
        tk.Entry(input_frame, textvariable=self.output_dir, width=50, bg="white").grid(row=1, column=1, padx=5)

        # --- LOG BOX (The "Matrix" Look) ---
        self.log_box = tk.Text(root, height=12, bg="black", fg="#00ff00", font=("Courier", 9))
        self.log_box.pack(padx=10, pady=5, fill="both")

        # --- BUTTONS ---
        btn_f = tk.Frame(root, bg="#d9d9d9")
        btn_f.pack(pady=10)

        tk.Button(btn_f, text="[ RUN_TRANSLATE ]", command=self.do_translate, width=20, bd=4).pack(side="left", padx=5)
        tk.Button(btn_f, text="[ RUN_DIGEST ]", command=self.do_digest, width=20, bd=4).pack(side="left", padx=5)
        tk.Button(btn_f, text="[ CLEAR_LOG ]", command=lambda: self.log_box.delete('1.0', tk.END), width=15, bd=4).pack(side="left", padx=5)

    def log(self, msg):
        self.log_box.insert(tk.END, f"> {msg}\n")
        self.log_box.see(tk.END)
        self.root.update()

    def select_file(self):
        path = filedialog.askopenfilename()
        if path: self.fasta_path.set(path)

    def select_folder(self):
        path = filedialog.askdirectory()
        if path: self.output_dir.set(path)

    def do_translate(self):
        input_f = self.fasta_path.get()
        out_root = self.output_dir.get()
        
        if not input_f or not out_root:
            messagebox.showerror("ERROR", "PATHS NOT DEFINED")
            return

        base = os.path.basename(input_f).split('.')[0]
        organism_folder = os.path.join(out_root, f"{base}_full_proteins")
        os.makedirs(organism_folder, exist_ok=True)

        self.log(f"INITIATING TRANSLATION: {base}")
        
        try:
            count = 0
            for record in SeqIO.parse(input_f, "fasta"):
                # DNA -> Protein string conversion
                full_protein = record.seq.translate()
                all_proteins = str(full_protein).split('*')
                real_proteins = [p for p in all_proteins if len(p) > 20]

                for i, p_seq in enumerate(real_proteins):
                    try:
                        analysed_p = ProteinAnalysis(p_seq)
                        weight = analysed_p.molecular_weight()
                        
                        file_name = f"{base}_fragment_{i}.txt"
                        file_path = os.path.join(organism_folder, file_name)

                        with open(file_path, 'w') as f:
                            f.write(f"Source: {record.description}\n")
                            f.write(f"Weight: {weight:.2f} Daltons\n")
                            f.write("-" * 20 + "\n")
                            f.write(p_seq)
                        count += 1
                    except: continue 
            
            self.log(f"SUCCESS: {count} PROTEINS IN {organism_folder}")
            messagebox.showinfo("SUCCESS", "TRANSLATION FINISHED")
        except Exception as e:
            self.log(f"CRITICAL ERROR: {str(e)}")

    def do_digest(self):
        input_f = self.fasta_path.get()
        out_root = self.output_dir.get()
        
        if not input_f or not out_root:
            messagebox.showerror("ERROR", "PATHS NOT DEFINED")
            return

        base = os.path.basename(input_f).split('.')[0]
        final_dir = os.path.join(out_root, f"{base}_digest_analysis")
        os.makedirs(final_dir, exist_ok=True)

        self.log(f"INITIATING DIGEST: {base}")
        
        try:
            count = 0
            for record in SeqIO.parse(input_f, "fasta"):
                # to_stop=False ensures we read the whole DNA sequence even past the first stop
                protein_string = str(record.seq.translate(to_stop=False))
                fragments = protein_string.split('*')
                real_frags = [f for f in fragments if len(f) > 0] # Catch all fragments
                
                self.log(f"SCANNING {record.id}...")

                for i, p_seq in enumerate(real_frags):
                    try:
                        analysed = ProteinAnalysis(p_seq)
                        weight = analysed.molecular_weight()
                        gravy = analysed.gravy()
                        location = "External" if gravy < 0 else "Internal"
                        
                        file_path = os.path.join(final_dir, f"part_{i}.txt")
                        with open(file_path, 'w') as f:
                            f.write(f"Source: {record.id}\n")
                            f.write(f"Weight: {weight:.2f} Da\n")
                            f.write(f"GRAVY: {gravy:.2f}\n")
                            f.write(f"Location: {location}\n")
                            f.write("-" * 20 + "\n")
                            f.write(p_seq)
                        count += 1
                    except: continue
            
            self.log(f"DIGEST COMPLETE: {count} FILES CREATED")
            messagebox.showinfo("SUCCESS", "DIGEST FINISHED")
        except Exception as e:
            self.log(f"CRITICAL ERROR: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = BioWorkbench(root)
    root.mainloop()