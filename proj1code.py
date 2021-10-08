#dna2rna reverse-complement rna2codon
def rna2dna(rna):
    dna = ''
    for symbol in rna:
        if symbol == 'A': #a to t
            dna = dna + 'T'
        elif symbol == 'U': #u to a
            dna = dna + 'A'
        elif symbol == 'C': #c to g
            dna = dna + 'G'
        elif symbol == 'G': #g to c
            dna = dna + 'C'
    return dna

def rna2codon(triplet):
    genetic_code = { #dict for codons.
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    allowed_codons = set('ACGU')
    if genetic_code.get(triplet) == None: #if it's not a valid triplet, the triplet is invalid.
        return "Invalid"
    return genetic_code.get(triplet)

def reverse_complement(dna):
    reverse = dna[::-1] #reverses string.
    revcomp = '' #blank reverse complement string.
    for i in reverse:
        if i == 'A': #a->t
            revcomp += 'T'
        elif i == 'T': #t->a
            revcomp += 'A'
        elif i == 'G': #g->c
            revcomp += 'C'
        elif i == 'C': #c->g
            revcomp += 'G'
        else:
            # returns an error if there is an erroneous entry in the dna string.
            return 'error in dna string: non-nucleotide detected.'
    return revcomp
