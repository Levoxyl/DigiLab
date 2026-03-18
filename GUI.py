import tkinter as tk
from tkinter import filedialog, messagebox
import os

# Ensure the folder name and file names match your directory exactly (case-sensitive)
from Processing import Translation
from Processing import Digest

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Entrez

class BioWorkbench:
    def __init__(self, root):
        self.root = root
        self.root.title("GENOMIC_ANALYZER_V1.1")
        self.root.geometry("700x600")
        self.root.configure(bg="#d9d9d9")

        self.fasta_path = tk.StringVar()
        self.output_dir = tk.StringVar()

        # --- HEADER ---
        tk.Label(root, text="LAB_TERMINAL: DNA_DIGEST_SYSTEM", bg="#000080", fg="white", 
                 font=("Courier", 12, "bold"), relief="raised", bd=3).pack(fill="x", pady=5)

        # --- INPUTS ---
        input_frame = tk.Frame(root, bg="#d9d9d9", bd=2, relief="groove")
        input_frame.pack(fill="x", padx=10, pady=10)

        tk.Button(input_frame, text="FILE...", command=self.select_file, width=10).grid(row=0, column=0, padx=5, pady=5)
        tk.Entry(input_frame, textvariable=self.fasta_path, width=60, bg="white").grid(row=0, column=1, padx=5)

        tk.Button(input_frame, text="DIR...", command=self.select_folder, width=10).grid(row=1, column=0, padx=5, pady=5)
        tk.Entry(input_frame, textvariable=self.output_dir, width=60, bg="white").grid(row=1, column=1, padx=5)

        # --- LOG BOX ---
        self.log_box = tk.Text(root, height=15, bg="black", fg="#00ff00", font=("Courier", 9))
        self.log_box.pack(padx=10, pady=5, fill="both")

        # --- BUTTONS ---
        btn_f = tk.Frame(root, bg="#d9d9d9")
        btn_f.pack(pady=10)

        tk.Button(btn_f, text="[ RUN_TRANSLATE ]", command=self.do_translate, width=18, bd=4).pack(side="left", padx=2)
        tk.Button(btn_f, text="[ RUN_DIGEST ]", command=self.do_digest, width=18, bd=4).pack(side="left", padx=2)
        tk.Button(btn_f, text="[ IDENTIFY_VIRUS ]", command=self.do_db_search, width=18, bd=4).pack(side="left", padx=2)
        tk.Button(btn_f, text="[ CLEAR ]", command=lambda: self.log_box.delete('1.0', tk.END), width=10, bd=4).pack(side="left", padx=2)
        tk.Button(btn_f, text="[ ANALYZE_CHEMISTRY ]", command=self.analyze_chemistry, width=20, bd=4).pack(side="left", padx=5)

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

    def web_identify(self, sequence_id):
        Entrez.email = "your_email@example.com"  # <-- MANDATORY: Put your email here
        self.log(f"FETCHING INFO FOR: {sequence_id}")
        
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="gb", retmode="text")
            record = handle.read()
            handle.close()
            
            for line in record.split('\n'):
                if "DEFINITION" in line:
                    info = line.replace("DEFINITION", "").strip()
                    self.log(f"OFFICIAL NAME: {info}")
                    break
        except Exception as e:
            self.log(f"API ERROR: {str(e)}")

    def do_db_search(self):
        input_f = self.fasta_path.get()
        if not input_f:
            messagebox.showerror("ERROR", "SELECT FASTA FILE FIRST")
            return
        
        try:
            # "Sniff" the ID out of the file automatically
            for record in SeqIO.parse(input_f, "fasta"):
                self.log(f"ID DETECTED: {record.id}")
                self.web_identify(record.id)
                break # Only need the first ID for identification
        except Exception as e:
            self.log(f"READ ERROR: {str(e)}")
            
    def analyze_chemistry(self):
        input_f = self.fasta_path.get()
        if not input_f:
            messagebox.showerror("ERROR", "SELECT FASTA FIRST")
            return
        
        self.log("RUNNING BIO-CHEMICAL SCAN...")
        
        try:
            for record in SeqIO.parse(input_f, "fasta"):
                protein_string = str(record.seq.translate(to_stop=False))
                fragments = [f for f in protein_string.split('*') if len(f) > 30]
                
                sticky_boss = max(fragments, key=lambda p: ProteinAnalysis(p).gravy())
                
                analysed = ProteinAnalysis(sticky_boss)
                self.log(f"--- RESULTS FOR {record.id} ---")
                self.log(f"STICKIEST PROTEIN: {len(sticky_boss)} amino acids")
                self.log(f"GRAVY SCORE: {analysed.gravy():.2f}")
                
                if analysed.gravy() > 0:
                    self.log("PREDICTION: Internal / membrane-Anchored")
                else:
                    self.log("PREDICTION: External / Surface-Exposed")
                    
                aromatic = analysed.aromaticity()
                self.log(f"STABILITY (Aromaticity): {aromatic*100:.1f}%")
                
        except Exception as e:
            self.log(f"ANALYSES ERROR: {str(e)}")
            
    def get_amino_distro(self):
        input_f = self.fasta_path.get()
        if not input_f:
            messagebox.showerror("ERROR", "SELECT FASTA FIRST")
            return

        self.log("CALCULATING AMINO ACID DISTRIBUTION...")
        
        try:
            for record in SeqIO.parse(input_f, "fasta"):
                full_p = record.seq.translate(to_stop=True)
                analysed = ProteinAnalysis(str(full_p))
                distro = analysed.get_amino_acids_percent()
                top_3 = sorted(distro.items(), key=lambda x: x[1], reverse=True)[:3]
                
                self.log(f"--- {record.id} COMPOSITION ---")
                for aa, percent in top_3:
                    self.log(f"AA '{aa}': {percent*100:.2f}%")
                
                if distro.get('L', 0) > 0.10:
                    self.log("NOTICE: High Leucine content detected.")
        except Exception as e:
            self.log(f"DISTRO ERROR: {str(e)}")

    # Update your CLEAR button in __init__ to this for speed:
    # tk.Button(btn_f, text="[ CLEAR ]", command=self.clear_log, width=10, bd=4).pack(side="left", padx=2)

    def clear_log(self):
        self.log_box.delete('1.0', tk.END)
        self.root.update_idletasks() # Forces an immediate visual refresh

    def do_translate(self):
        input_f = self.fasta_path.get()
        out_root = self.output_dir.get()
        
        if not input_f or not out_root:
            messagebox.showerror("ERROR", "SELECT FILE AND OUTPUT DIR")
            return
            
        self.log(f"STARTING TRANSLATION -> {out_root}")
        Translation.process_lab_directory(input_f, out_root)
        self.log("TRANSLATION COMPLETE.")

    def do_digest(self):
        input_f = self.fasta_path.get()
        out_root = self.output_dir.get()
        
        if not input_f or not out_root:
            messagebox.showerror("ERROR", "SELECT FILE AND OUTPUT DIR")
            return
            
        self.log(f"STARTING DIGEST -> {out_root}")
        Digest.find_virus_parts(input_f, out_root)
        self.log("DIGEST COMPLETE.")

if __name__ == "__main__":
    root = tk.Tk()
    app = BioWorkbench(root)
    root.mainloop()