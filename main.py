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

    retrived_nucleotides = bio_data_vis.retrive_first_nucleotides( file_path=r'./retrived_seqs.fna', n_nucleotides=13556 )
    
    ls_retrived_data = []

    if ( bio_data_vis._DBG0_ ):
        for k,v in retrived_nucleotides.items():
            retrived_nucleotides[k] = str(v)

        print( f"length of sequences: {len( retrived_nucleotides['NZ_AP022841.1'] )}" )

    # objects
    main_analysis = bio_data_vis.Analysis( r'./aligned_sequences.fna', 'clustal' )
    main_tests = bio_data_vis.Testing()

    bio_test_data = bio_data_vis.fasta_msa( retrived_nucleotides, 13556 )    
    id_indels = main_analysis.id_indels()
    aligned_seqs = bio_data_vis.create_custom_msa_data( bio_test_data, 13556, 'test_msa_file.fna' )

    main_analysis.vis_heatmap_msa( aligned_seqs )           # after MSA
    main_analysis.vis_heatmap_test( retrived_nucleotides )  # before MSA
    
    # indel_df = pd.DataFrame.from_dict( id_indels, orient='index', columns=['Gap count'] )
    
    # main_analysis.generate_phylo_tree( 'to_align.dnd', 'newick' )

    # print( f'retrived data:\n{ retrived_nucleotides[1] }' )
    ls_nucleotide_diversity = main_analysis.calc_nucleotide_diversity( retrived_nucleotides, 13556 )
    main_analysis.vis_genetic_diversity( ls_nucleotide_diversity )
    
    # print( f'nucleotide_diversity: { nucleotide_diversity }' )