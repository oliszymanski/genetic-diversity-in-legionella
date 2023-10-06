import data.bio_data_vis as bio_data_vis

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
    
    bio_data_vis.write_sequences_to_align( main_file_path, wanted_ids )
    

    # bio_data_vis.align_two_seq( 'global', ls_record_penumophila[0], ls_record_antarctica[0] )