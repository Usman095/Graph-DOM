import pandas as pd

from src import config, pathways, outputs, utils, families

if __name__ == '__main__':
    df = pd.read_excel(config.input_file_path)
    is_precursor = abs(df['Precursor m/z'] - df['fragments m/z']) < config.nominal_tolerance
    precursor_data = df[is_precursor]

    row_list = list(precursor_data['Chemical formula'])
    row_list.sort(key=lambda x: utils.get_mass(utils.get_formula(x)))
    row = [utils.get_formula(x) for x in row_list]
    masses = [utils.get_mass(x) for x in row]
    print(len(df))
    print(len(precursor_data))

    nominal_dict = {}
    for _, l_row in precursor_data.iterrows():
        if l_row['Chemical formula'] not in nominal_dict:
            nominal_dict[l_row['Chemical formula']] = l_row['Precursor m/z']


    df = df.sort_values(['Precursor m/z', 'fragments m/z'], ascending=[True, False]).reset_index(drop=True)
    groups = df.groupby('Precursor m/z', sort=False)

    pathway_dict, precursor_pathway_hist, o_count = pathways.generate_pathways_par(groups)

    print("Total number of precursors: {}".format(len(pathway_dict)))
    outputs.write_pathway_to_csv(pathway_dict)

    outputs.vk(row)
    outputs.core_dist_over_precursor(pathway_dict, x_axis="pre_id")
    outputs.pathway_dist_over_precursor(pathway_dict, precursor_pathway_hist, x_axis="pre_id")
    outputs.pathway_dist_over_oxygen_class(o_count)
    outputs.core_dist_over_oxygen_class(pathway_dict)

    pathway_dict_df = pd.DataFrame([[pre, core, index, set(row["path"].copy()), row["path"].copy(), row["CoreMass"]] 
                            for pre in pathway_dict.keys() 
                            for core in pathway_dict[pre].keys()
                            for index, row in pathway_dict[pre][core].items()],
                            columns=["Precursor", "Core-Fragment", "ID", "Path-Set", "Pathway", "Core-Mass"])
    pathway_dict_df = pathway_dict_df.astype({'Core-Mass': 'int32'})
    pathway_dict_df.drop(columns=["Path-Set"], inplace=True, axis=1)
    pathway_dict_df["Pre-Mass"] = [utils.get_mass(utils.get_formula(precursor)) for precursor in pathway_dict_df["Precursor"]]
    pathway_dict_df.sort_values(by=["Pre-Mass"], inplace=True, ascending=True)
    pathway_dict_df.reset_index(drop=True, inplace=True)

    ######## Families Start Here ##########
    print("Creating families...")
    roots, path_forest = families.get_path_forest(pathway_dict_df)
    family_dict = families.combine_families(roots, path_forest)
    outputs.write_families_to_csv(family_dict)

    print("Found {} total families.".format(len(family_dict)))

    outputs.isomers_vs_family_id(family_dict)
    outputs.write_cytoscape_family_graph(family_dict)
    outputs.family_parents_vs_oxygen_class(family_dict)
    outputs.family_size_dist(family_dict)

    ########### Stats ##############
    print("Total Number of Precursors: {}".format(len(set(list(precursor_data["Chemical formula"])))))
    pre_set = set()
    for prec_id, family_group in family_dict.items():
        for pre in prec_id:
            pre_set.add(pre)
    print("Precursors in Families: {}".format(len(pre_set)))

    total_cores, cores_in_families = outputs.core_coverage(nominal_dict, pathway_dict, family_dict)
    print("Total Number of Core-Fragments Identified: {}".format(total_cores))
    print("Cores-Fragments Covered in Families: {}".format(cores_in_families))

    total_frags, family_frags = outputs.fragment_coverage(nominal_dict, pathway_dict, family_dict)
    print("Total Number of Fragments in the Input: {}".format(total_frags))
    print("Fragments Covered in Families: {}".format(family_frags))