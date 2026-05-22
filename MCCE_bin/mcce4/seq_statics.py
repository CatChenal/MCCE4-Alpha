#!/usr/bin/env python3

"""
Module: seq_statics.py

Calculate the theoretical pI (isoelectric point) and charge of a
protein sequence using the callable SeqStatics class.
Call examples:
    theoret_pI, theoret_crg = SeqStatics(seq1, aa_code_len=1)()
    statics3 = SeqStatics(seq2, aa_code_len=3)()

"""
from collections import defaultdict


PHOS_RES = ["SEP","s", "TPO", "t", "PTR", "y"]
IONIZABLE_RES = ["ASP", "GLU", "ARG", "HIS", "LYS", "CYS", "TYR", "SEC",
                 # phosphorylated res:
                 "SEP", "TPO", "PTR",
                 # TER res: keep last!
                 "NTR", "CTR"]
ACIDIC_RES = ["ASP", "GLU", "SEC",
              # ... phosporylated res ...
             "SEP", "TPO",  "PTR"]
BASIC_RES = ["ARG", "HIS", "LYS"]

res3_to_res1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYD": "C",
    "CYL": "C",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "MSE": "M",
    "PHE": "F",
    "PRO": "P",
    "PYL": "O",
    "SEC": "U",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # ... phosporylated res ...
    "SEP": "s",  # phosphoSer   
    "TPO": "t",  # phosphoThr
    "PTR": "y",  # phosphoTyr
    # ... ligands ...
    "PL9": "MQ8",
    }

PREC = 2
NITER = 13
# NTR, CTR pKas adjusted to get crg=0 on polyA seq
TER_PKAS = [10.23, 3.0]

# ref: https://chem.libretexts.org/Courses/University_of_Arkansas_Little_Rock/CHEM_4320_5320%3A_Biochemistry_1/01%3A_Amino_Acids/1.4%3A_Reactions_of_Amino_Acids/1.4.1_Acid-base_Chemistry_of_Amino_Acids
SOLUTION_PKAS = {
    "ASP": 3.9,
    "GLU": 4.3,
    "CYS": 8.3,
    "ARG": 12.5,
    "HIS": 6.0,
    "LYS": 10.5,
    "SEC": 5.5,    # Selenocysteine, SEC -> U
    "TYR": 10.1,
    #... phospho res pKa: PMC6875518
    "SEP": [1.2, 6.0],  # phosphoSer
    "TPO": [1.3, 6.3],  # phosphoThr
    "PTR": [1.0, 5.8],  # phosphoTyr
}


class SeqStatics:
    def __init__(self, sequence: str, aa_code_len:int=1, n_iter: int=NITER):
        """Initialize SeqStatics instance with at least its sequence.
        Arguments:
         - sequence (str): Amino acid sequence (protein)
         - aa_code_len (int, 1, [1,3]): Indicates whether the AAs are given in 
           1-letter codes (default) or 3.
         - n_iter (int, NITER): Iterations for binary search.
        """
        if not sequence:
            raise ValueError("Empty sequence!")
        if aa_code_len not in [1, 3]:
            raise ValueError("aa_code_len must be either 1 or 3!")

        self.seq = sequence
        self.n_iter = n_iter

        if aa_code_len == 1:
            # convert dicts to 1-letter code:
            self.soln_pkas = dict((res3_to_res1.get(res), SOLUTION_PKAS.get(res))
                                  for res in SOLUTION_PKAS)
            # last 2: TER res: always included in get_theoretical_charge
            self.ionizable = [res3_to_res1.get(res) for res in IONIZABLE_RES[:-2]]
            Acids = [res3_to_res1.get(res) for res in ACIDIC_RES]
            Bases = [res3_to_res1.get(res) for res in BASIC_RES]
        else:
            # remove '(' & ')' from seq: 
            # case when entity_poly.pdbx_seq_one_letter_code downloaded from rcsb:
            try:
                ix = self.seq.index("(")
                self.seq = self.seq.replace("(","").replace(")","")
            except ValueError:
                pass
                
            self.soln_pkas = SOLUTION_PKAS
            self.ionizable = IONIZABLE_RES[:-2]
            Acids = ACIDIC_RES
            Bases = BASIC_RES
    
        ioniz_counts = defaultdict(int)
        self.phos_res = []
        for i in range(0, len(self.seq), aa_code_len):
            res = self.seq[i:i+aa_code_len]
            if res in self.ionizable:
                ioniz_counts[res] += 1
                if res in PHOS_RES:
                    self.phos_res.append(res)
        if not ioniz_counts:
            print("No ionizable residues in sequence: pI will be that of N-TER & C-TER.")

        self.ioniz_counts = dict(ioniz_counts)
        self.pk_dict = dict((res, self.soln_pkas.get(res)) for res in self.ioniz_counts)
        self.acids = [res for res in Acids if res in self.ioniz_counts]
        self.bases = [res for res in Bases  if res in self.ioniz_counts]

        return

    def get_theoretical_charge(self, ph: float) -> float:
        """Return the net charge at one pH.
        """
        # Positive charges
        charge = 1 / (1 + 10**(ph - TER_PKAS[0]))
        for b in self.bases:
            charge += self.ioniz_counts[b] / (1 + 10**(ph - self.pk_dict[b]))

        # Negative charges
        charge -= 1 / (1 + 10**(TER_PKAS[1] - ph))
        for a in self.acids:
            if a in self.phos_res:
                charge -= self.ioniz_counts[a] / (1 + 10**(self.pk_dict[a][0] - ph)) 
                charge -= self.ioniz_counts[a] / (1 + 10**(self.pk_dict[a][1] - ph))
            else:
                charge -= self.ioniz_counts[a] / (1 + 10**(self.pk_dict[a] - ph))

        return charge

    def get_theoretical_pI_crg(self, n_iter:int = None) -> list:
        """Binary search for pH where charge is 0.
        Increase iterations with n_iter (default:20) for higher precision.
        """
        if not self.soln_pkas:
            return [None, None]
    
        if n_iter is None:
            niter = self.n_iter
        else:
            niter = n_iter
    
        low, high = 0.0, 14.0
        for _ in range(niter):
            mid = (low + high) / 2
            if self.get_theoretical_charge(mid) > 0:
                low = mid
            else:
                high = mid
        
        return [round(mid, PREC), round(self.get_theoretical_charge(7.0), PREC)]

    def __call__(self, n_iter:int = None) -> list:
        """Call method to enable single command results, e.g.:
           `statics1 = SeqStatics(seq1_4lzt, aa_code_len=1)()`
        """
        return self.get_theoretical_pI_crg(n_iter=n_iter)
