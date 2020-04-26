from assignments.clustal_parser.clustal_parser import parse_clustal


with open('clustal_parser/test_data/clustal_assigned') as f:
    msa = parse_clustal(f)
    