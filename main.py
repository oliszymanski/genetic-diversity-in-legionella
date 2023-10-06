import data.bio_data_vis as bio_data_vis

from Bio import SeqIO

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
    
    
    # for record in SeqIO.parse( main_file_path, 'fasta' ):
    #     print( f'sequence id: { record.id }' )



    # ls_record_antarctica = [ str( seq ) for seq in ls_record_penumophila ]


    # if ( bio_data_vis._DBG0_ ):
    #     print( "ls_record_penumophila =", ls_record_penumophila )

    #     print( 'length of genome: ', len( ls_record_penumophila[ 0 ] ) )
        



    # seq_00 = SeqIO.parse(  )
    # seq_01 = ''

    # bio_data_vis.align_two_seq( 'global', ls_record_penumophila[0], ls_record_antarctica[0] )