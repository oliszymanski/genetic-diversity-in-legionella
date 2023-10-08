#========================================================
#   IMPORTS
#========================================================

import os

import matplotlib.pyplot as plt

from Bio import SeqIO, AlignIO, pairwise2, Phylo
from Bio.pairwise2 import format_alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline

# for debugging and testing
_DBG0_ = True



#========================================================
#   FUNCTIONS
#========================================================

def write_sequences_to_align( sequences_file_path : str, wanted_ids, file_type='fasta' ):
    """
    :param sequences_file_path: path of file with sequences,
    :param wanted_ids: list of ids of sequence records;
    """

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
    """
    :param file_name: name of file with sequences for MSA,
    :param file_type: type of file;
    """

    aligned_file = 'aligned_sequences.fna'
    clustalw_cline = ClustalwCommandline( 'C:/ClustalW2/clustalw2.exe', infile=file_name, outfile=aligned_file )
    clustalw_cline()

    aligned = AlignIO.read( aligned_file, file_type )
    aligned_seq = MultipleSeqAlignment( aligned )

    return aligned_seq



def retrive_first_nucleotides( file_path, n_nucleotides : int, file_type='fasta' ):
    """
    :param file_path: path of file (or file name in same dir),
    :param n_nucleotides: number of first n nucleotides;
    """

    seq_dict = {}

    for record in SeqIO.parse( file_path, file_type ):
        record_id = record.id
        seq = record.seq[ :n_nucleotides ]
        seq_dict[ record_id ] = seq

    return seq_dict



def write_fasta_file( in_dict : dict, output_file : str ):
    """
    :param in_dict: input dictonary with fasta data,
    :param output_file: output file after operations;
    """

    edit_file = open( output_file, 'w' )

    for record_id, seq in in_dict.items():
        edit_file.write( f'>{ record_id }\n{ seq }\n' )

    return None



def alignment_length( file_path: str, format : str ):
    """
    :param file_path: path of file (or file name in same dir),
    :param format: format of file;
    """

    alignment = AlignIO.read( file_path , format )
    alignment_length = alignment.get_alignment_length()

    print( 'Alignment length =', alignment_length )

    return alignment_length



def fasta_msa( file_path : str, file='type' ):
    return None




class Analysis():
    def __init__( self, aligned_sequences_file : str, file_type ):
        """
        :param aligned_sequences_file: file with aligned sequences,
        :param file_type: type of input file;
        """

        self.aligned_sequences_file = aligned_sequences_file
        self.alignments = AlignIO.read( self.aligned_sequences_file, file_type )

        return None
    


    def generate_phylo_tree(self, read_dnd_file : str, file_type='newick'):
        """
        :param read_dnd_file: dnd file from MSA file;
        """

        tree = Phylo.read( read_dnd_file, file_type )

        fig, ax = plt.subplots( figsize=(10, 10) )
        Phylo.draw( tree, axes=ax, do_show=False )
        plt.show()


        return None
    


    def id_indels( self ):

        indel_data = {}

        for i in range( len( self.alignments[0] ) ):
            col = self.alignments[ :, i ]
            gap_count = col.count( '-' )

            if (_DBG0_): print( f'column: {i}: {col}, gap count: {gap_count}' )

            if ( gap_count > 0 ):
                indel_data[i] = gap_count

        return indel_data
    


    def calc_nucleotide_diversity( self ) -> int:

        seq_length = len( self.alignments[0].seq )
        nucleotide_freq = [0] * seq_length

        for record in self.alignments:
            for i, nucleotide in enumerate( record.seq ):
                if ( nucleotide != '-' ): nucleotide_freq[i] += 1

        pi = 1 - sum( ( freq / len( self.alignments ) ) ** 2 for freq in nucleotide_freq )

        return pi