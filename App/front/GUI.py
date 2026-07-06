import tkinter as tk
from tkinter import filedialog, messagebox

import os
import sys

root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if root_path not in sys.path:
    sys.path.insert(0, root_path)

from pathlib import Path
from dotenv import load_dotenv  

from App.back.Processing import Translation
from App.back.Processing import Digest

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Entrez

env_path = Path(__file__).resolve().parent.parent.parent / '.env'
load_dotenv(dotenv_path=env_path)


class BioWorkbench:
    def __init__(self, root, theme_pack):
        self.root = root
        self.theme_pack = theme_pack

        self.theme = self.theme_pack.config
        self.btn_style = self.theme_pack.button

        self.fasta_path = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.log_box = None 

        self.ui = self.theme_pack.layout_class(self.root, self.theme_pack, controller=self)
        self.ui.build_interface()

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

    def clear_log(self):
        self.log_box.delete('1.0', tk.END)
        self.root.update_idletasks()

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