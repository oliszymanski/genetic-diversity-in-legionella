import data.bio_data_vis as bio_data_vis

from Bio import SeqIO, Phylo

if ( __name__ == '__main__' ):
    main_file_path = r'./sequences.fna'

    wanted_ids = [
        'NZ_AP022841.1',        # Legionella antarctica
        'NZ_CP015943.1',        # Legionella pneumophila
        'NZ_UGOW01000005.1',    # Legionella quateirensis
        'NC_003197.2',          # Salmonella enterica
        'NZ_LNYW01000001.1',    # Legionella shakespearei
        'NZ_LT906442.1'         # Legionella waltersii
    ]

    # retrived_nucleotides = bio_data_vis.retrive_first_nucleotides( file_path=r'./retrived_seqs.fna', n_nucleotides=60000 )
    # bio_data_vis.write_fasta_file( retrived_nucleotides, r'to_align.fna' )

    # aligned_sequences = bio_data_vis.align_multiple_seq_from_file( r'./to_align.fna' )

    main_analysis = bio_data_vis.Analysis( r'./aligned_sequences.fna', 'clustal' )
    id_indels = main_analysis.id_indels()
    main_analysis.generate_phylo_tree( 'to_align.dnd', 'newick' )

