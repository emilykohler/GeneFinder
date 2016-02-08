# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE - wat 

@author: Emily Kohler

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'

    The unit tests check whether the function is working properly and producing the 
    complementary nucleotide.

    """
    X = nucleotide
    if X == 'A':
        return 'T'
    elif X == 'T':
        return 'A'
    elif X == 'C':
        return 'G'
    elif X == 'G':
        return 'C'
    else:
        return 'Invalid Input'

    # TODO: implement this
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'

    Unit Tests are used to check whether the function is working, by providing strings of 
    different lengths and testing for the reverse complement (where the string is first 
    reversed and the the complements are found).

    """
    string = ''
    length = len(dna)
    reverse = dna[::-1]
    for value in reverse:
        complement = get_complement(value)
        string = string + complement
    return string

    # TODO: implement this
    pass

#print(get_reverse_complement('ATGCCCGCTTT'))


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF('ATGAAT')
    'ATGAAT'

    Tests for the function to cut of the sequence at the right time and whehter the function
    reuturns the original string if a stop codon isn't found.

    """
# stop codons are: tga, tag, taa
    length = len(dna)
    end_index = 0
    start_index = -3
    ORF = ''
    while end_index <= length:
        end_index = end_index + 3
        start_index = start_index + 3
        new_string = dna[start_index:end_index]
        if new_string == 'TAG' or new_string == 'TGA' or new_string == 'TAA':
            return ORF
        else:
            ORF = ORF + new_string
    return dna
        



    # TODO: implement this
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']

    >>> find_all_ORFs_oneframe("CATGAATGTAGATAGATGTGCCC")
    ['ATGTGCCC']

    Checks to see whether the find_all_ORFs_oneframe function first searches for an "ATG" and then uses 
    the rest_of_ORF function to find the ORFs.

    """
    length = len(dna)
    end_index = 0
    start_index = -3
    ORF = []
    indices = []
    indices2 = []
    while end_index <= length:
        end_index = end_index + 3
        start_index = start_index + 3
        new_string = dna[start_index:end_index]
        if new_string == 'ATG':
            indices.append(start_index)
        elif new_string == 'TAG' or new_string == 'TGA' or new_string == 'TAA' :
            indices2.append(start_index)
    #indices contain where the start codons are
    #print indices
    #print indices2

    for stop in indices2:
        for start in indices[indices2.index(stop)+1:]:
            if start < stop:
                indices.remove(start)

    for value in indices:
        string1 = dna[value:]
        string2 = rest_of_ORF(string1)
        ORF.append(string2)
    return ORF





    # TODO: implement this
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs('AAATTT')
    []

    If the dna sequence does not contain an 'ATG', there are no ORFs and the function returns
    an empty list.

    """
    value = -1
    indices = []
    indices2 = []
    while value <= 1:
        value = value + 1
        dna1 = dna[value:]
        indices.append(dna1)
    for string in indices:
        string1 = find_all_ORFs_oneframe(string)
        indices2 += string1
    return indices2


    # TODO: implement this
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGATAGCATCAAA")
    ['ATGCGA', 'ATGCTATCGCAT']
    >>> find_all_ORFs_both_strands("AAAAAA")
    []

    Finds the ORFs in the dna strand and then the reverse dna strand. Also, when there is 
    not an "ATG" the functions returns an empty list.
    """

    dna_reverse_complement = get_reverse_complement(dna)
    ORFs1 = find_all_ORFs(dna)
    ORFs2 = find_all_ORFs(dna_reverse_complement)
    ORFs3 = ORFs1 + ORFs2
    return ORFs3
   
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    dna1 = find_all_ORFs_both_strands(dna)
    length = len(dna1)
    max_length = 0
    max_strand = ''
    if length == 1:
        return str(dna1[0])
    else:
        for i in dna1:
            if len(i) >= max_length:
                max_length = len(i)
                max_strand = i
    return max_strand


    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """

    shuffled_dna_list = []
    ORF_length_list = []
    i  = 0
    while i <= num_trials:
        i = i + 1
        shuffled_dna = shuffle_string(dna)
        shuffled_dna_list.append(shuffled_dna)
    for dna_strand in shuffled_dna_list:
        non_coding_ORF = longest_ORF(dna_strand)
        ORF_length_list.append(len(non_coding_ORF))
    return max(ORF_length_list)

    
    # TODO: implement this
    pass
    


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    length = len(dna)
    end_index = 0
    start_index = -3
    amino_acid_string = ''
    while end_index <= length - 3:
        end_index = end_index + 3
        start_index = start_index + 3
        new_string = dna[start_index:end_index]
        amino = aa_table[new_string]
        amino_acid_string = amino_acid_string + amino
    return amino_acid_string


    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = longest_ORF_noncoding(dna, 1500)
    all_dna = find_all_ORFs_both_strands(dna)
    amino_list = []
    for dna_sequence in all_dna:
            if len(dna_sequence) >= threshold:
                amino_sequence = coding_strand_to_AA(dna_sequence)
                amino_list.append(amino_sequence)
    return amino_list


    # TODO: implement this
    pass

#if __name__ == "__main__":
   # import doctest
   # doctest.testmod()

#print(get_reverse_complement('ATGCCCGCTTT'))

#print(rest_of_ORF('ATGCATGAATGTAGATAGATGTGCCC'))
#print(rest_of_ORF('ATGAAT'))

#print(find_all_ORFs_oneframe("ATGCATATGGAATGTAGATAGATGTGCACCATGTAG"))

#print(find_all_ORFs("ATGCATGAATGTAG"))

#print(find_all_ORFs_both_strands("AAAAAA"))

#print(longest_ORF("ATGCGAATGTAGCATCAAA"))


#print(longest_ORF_noncoding('ATGCGAATGTACCCGTAGTAGCATCAAA', 6))

#print(coding_strand_to_AA("ATGCCCGCTTT"))

from load import load_seq
dna = load_seq("./data/X73525.fa")
print(gene_finder(dna))