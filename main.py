import data.bio_data_vis as bio_data_vis

if (__name__ == '__main__'):
    # data_vis.align_two_seq( type='local', alignment_00="ACCGT", alignment_01="ACG" )

    aligned_seq_records = bio_data_vis.align_multiple_seq( './data/genomes/', format='fasta' )

    
    # print( seq[0] ) # sequence of the first file .fna file