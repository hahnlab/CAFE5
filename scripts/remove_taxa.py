import sys

def get_taxa_cols(header_order_dict, taxa_list):
    taxa_cols = [header_order_dict[taxon_name] for taxon_name in taxa_list]
    return taxa_cols
    
def remove_taxa(input_file_path, output_file_path, taxa_to_remove):
    taxa_list = taxa_to_remove.split(',')
    header_order = dict() # node name:idx
    node_name_idx = 0

    with open(output_file_path, 'w') as output_file:
        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                line = line.rstrip()

                if line.startswith('#'):
                    node_name = line.lstrip('#')

                    if not node_name in taxa_list:
                        output_file.write(line + '\n')

                    header_order[node_name] = node_name_idx
                    node_name_idx += 1

                else:
                    taxa_cols_remove = get_taxa_cols(header_order, taxa_list)
                    tokens = line.split('\t')
                    tokens_to_print = [token for idx, token in enumerate(tokens) if not idx in taxa_cols_remove]
                    output_file.write('\t'.join(tokens_to_print) + '\n')
                    

if __name__ == "__main__":
    # sys.argv[1]: cafexp input file
    # sys.argv[2]: output file name
    # sys.argv[3]: taxa names, separated by , (to remove)
    
    remove_taxa(sys.argv[1], sys.argv[2], sys.argv[3])
