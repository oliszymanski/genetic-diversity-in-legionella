#========================================================
#   IMPORTS
#========================================================

import os

# import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline

# for debugging and testing
_DBG0_ = False



#========================================================
#   FUNCTIONS
#========================================================

def write_sequences_to_align( sequences_file_path : str, wanted_ids, file_type='fasta' ):

    wanted_records = []

    records = SeqIO.to_dict( SeqIO.parse( sequences_file_path, file_type ) )

    for record_id in wanted_ids:
        if ( record_id in records ):
            wanted_records.append( records[record_id] )

    SeqIO.write( wanted_records, 'retrived_seqs.fna', file_type )


    print( f"wanted seqs:\n{ wanted_records }" )


    return None



def align_two_seq( type : str, seq_00, seq_01 ):
    """
    :param type: define type of alignment'
    :param alignment_00: first sequence to align'
    :param alignment_01: second sequence to align'
    :return: formatted_alignment;
    """

    if ( type == 'global' ):
        alignments = pairwise2.align.globalxx( seq_00, seq_01 )
        formatted_alignment = format_alignment( *alignments[0] )

        if (_DBG0_): print( formatted_alignment )

        return formatted_alignment


    elif ( type == 'local' ):
        alignments = pairwise2.align.localxx( seq_00, seq_01 )
        formatted_alignment = format_alignment( *alignments[0] )

        if (_DBG0_): print( formatted_alignment )

        return formatted_alignment



def align_multiple_seq_from_dir( genomes_dir : str, format: str ):
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

    aligned = AlignIO.read( aligned_file, 'fasta' )
    aligned = MultipleSeqAlignment( aligned )

    return aligned



def align_multiple_seq_from_file( file_name : str, file_type='fasta' ):

    aligned_file = 'aligned_sequences.fna'
    clustalw_cline = ClustalwCommandline( 'C:/ClustalW2/clustalw2.exe', infile=file_name, outfile=aligned_file )
    clustalw_cline()

    aligned = AlignIO.read( aligned_file, file_type )
    aligned_seq = MultipleSeqAlignment( aligned )

    return aligned_seq



def alignment_length( file_path: str, format : str ):

    alignment = AlignIO.read( file_path , format )
    alignment_length = alignment.get_alignment_length()

    print( 'Alignment length =', alignment_length )

    return alignment_length


