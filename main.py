#========================================================
#   IMPORTS
#========================================================

import data.bio_data_vis as bio_data_vis



#========================================================
#   MAIN
#========================================================

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

    retrived_nucleotides = bio_data_vis.retrive_first_nucleotides( file_path=r'./retrived_seqs.fna', n_nucleotides=60000 )
    
    if ( bio_data_vis._DBG1_ ):
        for k,v in retrived_nucleotides.items():
            print( v[1412]  )
    
    # bio_data_vis.write_fasta_file( retrived_nucleotides, r'to_align.fna' )

    # aligned_sequences = bio_data_vis.align_multiple_seq_from_file( r'./to_align.fna' )

    # objects
    main_analysis = bio_data_vis.Analysis( r'./aligned_sequences.fna', 'clustal' )
    main_tests = bio_data_vis.Testing()

    bio_test_data = bio_data_vis.fasta_msa( retrived_nucleotides, 13556 )
    
    main_tests.test_alignments( retrived_nucleotides, bio_test_data )
    
    id_indels = main_analysis.id_indels()

    

    
    # indel_df = pd.DataFrame.from_dict( id_indels, orient='index', columns=['Gap count'] )
    
    # main_analysis.generate_phylo_tree( 'to_align.dnd', 'newick' )

    # nucleotide_diversity = main_analysis.calc_nucleotide_diversity()
    # print( f'nucleotide_diversity: { nucleotide_diversity }' )