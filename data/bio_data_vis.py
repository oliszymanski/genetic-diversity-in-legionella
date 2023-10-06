#========================================================
#   IMPORTS
#========================================================

import os

# import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, pairwise2, Phylo
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# for debugging and testing
_DBG0_ = True



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



def retrive_first_nucleotides( file_path, n_nucleotides : int, file_type='fasta' ):
    
    seq_dict = {}

    for record in SeqIO.parse( file_path, file_type ):
        record_id = record.id
        seq = record.seq[ :n_nucleotides ]
        seq_dict[ record_id ] = seq

    return seq_dict



def write_fasta_file( in_dict : dict, output_file : str ):

    edit_file = open( output_file, 'w' )

    for record_id, seq in in_dict.items():
        edit_file.write( f'>{ record_id }\n{ seq }\n' )

    return None



def alignment_length( file_path: str, format : str ):

    alignment = AlignIO.read( file_path , format )
    alignment_length = alignment.get_alignment_length()

    print( 'Alignment length =', alignment_length )

    return alignment_length



class Analysis():
    def __init__( self, aligned_sequences_file : str, file_type ):
        
        self.aligned_sequences_file = aligned_sequences_file
        self.file_type = file_type
        self.alignments = AlignIO.read( self.aligned_sequences_file, file_type )

        return None
    


    def generate_phylo_tree(self, output_file_name : str, output_file_type : str):
        calc = DistanceCalculator( 'identity' )
        distance = calc.get_distance( self.alignments )

        constructor = DistanceTreeConstructor( calc, method='nj' )
        tree = constructor.build_tree( distance )

        Phylo.write( tree, output_file_name, output_file_type )


        return None
    

    def id_indels( self ):

        # indel_positions = []
        # indel_lengths = []

        indel_data = {}

        for i in range( len( self.alignments[0] ) ):
            col = self.alignments[ :, i ]
            gap_count = col.count( '-' )

            if ( gap_count > 0 ):
                indel_data[i] = gap_count

        return indel_data