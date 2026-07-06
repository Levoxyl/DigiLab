from .Digest import Digest
from .Translation import Translation

class ProcessManager:
    """Acts as a coordinator master class between the processing units of the backend."""

    @staticmethod
    def parse_dna(fasta_path):
        print(f"[PROCESSING ENGINE] Parsing sequence via manager from: {fasta_path}")
        return "Parsed structural data from submodules"