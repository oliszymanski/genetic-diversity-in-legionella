#========================================================
#   IMPORTS
#========================================================

import os

import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
# from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalwCommandline

# for debugging and testing
_DBG0_ = True



#========================================================
#   FUNCTIONS
#========================================================

def say_hi():
    print('hello there')

    return



def align_two_seq( type : str, alignment_00, alignment_01 ):
    """
    :param type: define type of alignment'
    :param alignment_00: first sequence to align'
    :param alignment_01: second sequence to align'
    :return: formatted_alignment;
    """

    if ( type == 'global' ):
        alignments = pairwise2.align.globalxx( alignment_00, alignment_01 )
        formatted_alignment = format_alignment( *alignments[0] )

        if (_DBG0_): print( formatted_alignment )

        return formatted_alignment


    elif ( type == 'local' ):
        alignments = pairwise2.align.localxx( alignment_00, alignment_01 )
        formatted_alignment = format_alignment( *alignments[0] )

        if (_DBG0_): print( formatted_alignment )

        return formatted_alignment



def align_multiple_seq( genomes_dir : str, format: str ):
    """
    :param genomes_dir: directory containing all genomes to analyse,
    :param format: file format (set by default: fasta file),
    :return: multiple alignment;
    """

    seq_records = []
    temp_file = 'sequences.fna'
    
    for filename in os.listdir( genomes_dir ):
        if ( filename.endswith( '.fna' ) or filename.endswith( '.fasta' ) ):    # by default .fna and .fasta (change in the future)
            single_file_path = os.path.join( genomes_dir, filename )            # generates path to a single file

            for record in SeqIO.parse( single_file_path, format ):
                seq_records.append( record )


    to_align = 'sequences.fna'
    aligned_file = 'aligned_sequences.fna'

    SeqIO.write( seq_records, to_align, 'fasta' )


    clustalw_cline = ClustalwCommandline( 'C:/ClustalW2/clustalw2.exe', infile=to_align, outfile=aligned_file )
    clustalw_cline()

    # aligned = AlignIO.read( aligned_file, 'fasta' )

    return seq_records




def alignment_length( file_path: str, format : str ):

    alignment = AlignIO.read( file_path , format )
    alignment_length = alignment.get_alignment_length()

    print( 'Alignment length =', alignment_length )

    return alignment_length