import csv
from collections import defaultdict
from os.path import join

import matplotlib
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from src import utils

plot_path = 'output/plots'
file_path = 'output/files'


def vk(mol_dicts):
    hoc = []
    ooc = []
    for dic in mol_dicts:
        if 0.2 < dic['O'] / dic['C'] < 1.0:
            hoc.append(dic['H'] / dic['C'])
            ooc.append(dic['O'] / dic['C'])
    
    plt.scatter(ooc, hoc, s=5, c='green')
    plt.ylabel('H/C', fontsize=20)
    plt.xlabel('O/C', fontsize=20)
    plt.title('O/C H/C plot')
    plt.savefig(join(plot_path, 'Van_Krevelen_all_classes.png'), dpi=600)
    plt.close()


def core_dist_over_precursor(pathway_dict, x_axis='pre_mz'):
    pre_mass_num_core = []
    num_cores = []
    for pre, core_dict in pathway_dict.items():
        pre_mass_num_core.append([utils.get_mass(utils.get_formula(pre))/1000, len(core_dict), pre])
        
    pre_mass_num_core.sort(key=lambda x: x[0])
    pre_mz, num_cores, pre = map(list, zip(*pre_mass_num_core))

    print("Number of Core-Fragments in Pathways: {}".format(sum(num_cores)))
    if x_axis == 'pre_mz':
        plt.bar(pre_mz, num_cores,color='green')
        plt.xlabel('Precursor m/z')
    else:
        plt.bar(range(len(pre_mz)), num_cores,color='green')
        plt.xlabel('Precursor Id')
    plt.ylabel('Number of Core-Fragments')
    plt.title('Distribution of Core-Fragments over Precursors')
    plt.savefig(join(plot_path, 'Cores_vs_Precursors.png'), dpi=600)
    plt.close()
    f = open(join(file_path, 'Cores_vs_Precursors.csv'), "a")
    f.write("Precuresor ID, Precursor m/z, Number of Cores, Precursor Formula\n")
    for pre_id, (l_pre_mz, l_num_cores, l_pre) in enumerate(zip(pre_mz, num_cores, pre)):
        f.write("{}, {}, {}, {}\n".format(pre_id, l_pre_mz, l_num_cores, l_pre))
    f.close()


def pathway_dist_over_precursor(pathway_dict, precursor_pathway_hist, x_axis='pre_mz'):
    pre_mass_num_core = []
    num_cores = []
    for pre, core_dict in pathway_dict.items():
        pre_mass_num_core.append([utils.get_mass(utils.get_formula(pre))/1000, len(core_dict), pre])
        
    pre_mass_num_core.sort(key=lambda x: x[0])
    pre_mz, num_cores, pre = map(list, zip(*pre_mass_num_core))

    precursor_pathway_hist.sort(key = lambda x: x[0], reverse=False)
    x, y, z = map(list, zip(*precursor_pathway_hist)) # pre_mz, pre_formula, num_pathways

    print("Total Number of Pathways Identified: {}".format(sum(z)))    
    if x_axis == "pre_mz":
        plt.bar(x, z, color='green')
        plt.xlabel('Precursor m/z')
    else:
        plt.bar(range(len(x)), z, color='green')
        plt.xlabel('Precursor Id')
    plt.yscale("log")
    plt.ylabel('Number of Pathways')
    plt.suptitle('Pathway Count Distribution over Precursor')
    pre_path = join(plot_path, 'Pathway_Dist.png')
    plt.savefig(pre_path, dpi=600)
    plt.close()
    print('Saved plot at: "{}"'.format(pre_path))

    pre_path = join(file_path, 'Pathway_Dist.csv')
    f = open(pre_path, "a")
    f.write("Precursor ID, Precursor m/z, Number of Pathways, Precursor Formula\n")
    for pre_id, (l_pre_mz, l_path, l_pre) in enumerate(zip(x, z, y)):
        f.write("{}, {}, {}, {}\n".format(pre_id, l_pre_mz, l_path, l_pre))
    f.close()
    print('Saved file at: "{}"'.format(pre_path))


def pathway_dist_over_oxygen_class(o_count):
    print("Pathways Distribution over Oxygen Class:\n{}".format(o_count))
    x_pos = range(3, 16)
    labels = range(3, 16)
    plt.bar(x_pos, o_count[3: 16], align='center', alpha=0.5,color='green')
    plt.xticks(x_pos, labels)
    plt.ylabel('Number of Pathways')
    plt.xlabel('Oxygen Class')
    plt.title('Distribution of Pathways over Oxygen Class')
    o_class_path = join(plot_path, 'Pathways_vs_Oxygen_Class.png')
    plt.savefig(o_class_path, dpi=600)
    plt.close()
    print('Saved plot at: "{}"'.format(o_class_path))


def core_dist_over_oxygen_class(pathway_dict):
    o_len = 20
    o_count = [0] * o_len
    for pre, core_dict in pathway_dict.items():
        if 'N' not in pre:
            o_count[utils.get_formula(pre)["O"]] += len(core_dict)

    x_pos = range(3, 16)
    print(len(x_pos))
    labels = range(3, 16)
    print("Core Fragment Distribution over Oxygen Class:\n{}".format(o_count))
    plt.bar(x_pos, o_count[3: 16], align='center', alpha=0.5, color='green')
    plt.xticks(x_pos, labels)
    plt.ylabel('Number of Core-Fragments')
    plt.xlabel('Oxygen Class')
    plt.title('Distribtion of Core-Fragments over Oxygen Class')
    path = join(plot_path, 'Core-Fragment_vs_Oxygen_Class.png')
    plt.savefig(path, dpi=600)
    plt.close()
    print('Saved plot at: "{}"'.format(path))


def write_pathway_to_csv(pathway_dict):
    pathway_path = join(file_path, 'Fragmentation_Pathways.csv')
    with open(pathway_path, "w", newline='') as f:
        fieldnames = ["Precursor", "Core-Fragment", "ID", "Pathway", "CoreMass"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for pre, cores in pathway_dict.items():
            row = {}
            row["Precursor"] = pre
            for core, pathways in cores.items():
                row["Core-Fragment"] = core
                for ID, pathway in pathways.items():
                    row["ID"] = str(ID)
                    row["Pathway"] = pathway["path"]
                    row["CoreMass"] = pathway["CoreMass"]
                    #row.update(pathway)
                    writer.writerow(row)
                    row["Precursor"] = ""
                    row["Core-Fragment"] = ""
    
    print('Wrote Pathways to: "{}"'.format(pathway_path))


def write_families_to_csv(family_dict):
    family_output = []
    for prec_id, family_group in family_dict.items():
        # family_output.append([family_group, ""])
        for prec in prec_id:
            family_output.append([prec, "", "", "", "", "", ""])
        for family in family_group:
            for row in family:
                family_output.append([""] + row)
            family_output.append(["", "", "", "", "", "", ""])
        family_output.append(["", "", "", "", "", "", ""])
        family_output.append(["", "", "", "", "", "", ""])

    out_df = pd.DataFrame(family_output, columns=['Family-ID', 'Precursor', 'Core-Fragment', 'Pathway', 'Neutral-Loss', 'Overlap-Path', 'Short-Key'])
    # out_df = pd.DataFrame(family_output, columns=['Count', 'Family-ID'])
    family_path = join(file_path, 'Families.csv')
    out_df.to_csv(family_path)
    print('Wrote Families to: "{}"'.format(family_path))


def write_families_to_csv_short(family_dict):
    family_output = []
    for prec_id, family_group in family_dict.items():
        # family_output.append([family_group, ""])
        line = []
        count = 0
        for _ in family_group:
            count += 1
        line.append(count)
        prec_str = []
        for prec in prec_id:
            prec_str.append(prec)
        line.append(', '.join(prec_str))
        family_output.append(line)
        #     for row in family:
        #         family_output.append([""] + row)
        #     family_output.append(["", "", "", "", "", "", ""])
        # family_output.append(["", "", "", "", "", "", ""])
        # family_output.append(["", "", "", "", "", "", ""])

    # out_df = pd.DataFrame(family_output, columns=['Family-ID', 'Precursor', 'Core-Fragment', 'Pathway', 'Neutral-Loss', 'Overlap-Path', 'Short-Key'])
    out_df = pd.DataFrame(family_output, columns=['# of Overlap Pathways', 'Family'])
    family_path = join(file_path, 'Families-Short.csv')
    out_df.to_csv(family_path)
    print('Wrote Families to: "{}"'.format(family_path))


def write_fam4_to_csv(family_dict):
    '''Write only families of size >= 4 to csv'''
    family_output = []
    for prec_id, family_group in family_dict.items():
        # family_output.append([family_group, ""])
        if len(prec_id) >= 4:
            for prec in prec_id:
                family_output.append([prec, "", "", "", "", "", ""])
            for family in family_group:
                for row in family:
                    family_output.append([""] + row)
                family_output.append(["", "", "", "", "", "", ""])
            family_output.append(["", "", "", "", "", "", ""])
            family_output.append(["", "", "", "", "", "", ""])

    out_df = pd.DataFrame(family_output, columns=['Family-ID', 'Precursor', 'Core-Fragment', 'Pathway', 'Neutral-Loss', 'Overlap-Path', 'Short-Key'])
    # out_df = pd.DataFrame(family_output, columns=['Count', 'Family-ID'])
    family_path = join(file_path, 'Families-size-4-or-more.csv')
    out_df.to_csv(family_path)
    print('Wrote Families to: "{}"'.format(family_path))


def write_fam5_to_csv(family_dict):
    '''Write only families of size >= 4 to csv'''
    family_output = []
    for prec_id, family_group in family_dict.items():
        # family_output.append([family_group, ""])
        if len(prec_id) >= 5:
            for prec in prec_id:
                family_output.append([prec, "", "", "", "", "", ""])
            for family in family_group:
                for row in family:
                    family_output.append([""] + row)
                family_output.append(["", "", "", "", "", "", ""])
            family_output.append(["", "", "", "", "", "", ""])
            family_output.append(["", "", "", "", "", "", ""])

    out_df = pd.DataFrame(family_output, columns=['Family-ID', 'Precursor', 'Core-Fragment', 'Pathway', 'Neutral-Loss', 'Overlap-Path', 'Short-Key'])
    # out_df = pd.DataFrame(family_output, columns=['Count', 'Family-ID'])
    family_path = join(file_path, 'Families-size-4-or-more.csv')
    out_df.to_csv(family_path)
    print('Wrote Families to: "{}"'.format(family_path))


def isomers_vs_family_id(family_dict):
    dist = []
    for prec_id, family_group in family_dict.items():
        dist.append(len(family_group))
    #     if len(family_group) > 20000:
    #         print(prec_id)
    dist.sort(reverse=True)
    print("Isomers vs family id:")
    print(dist)
    # dist = dist[:20]
    plt.bar(range(len(dist)), dist, align='center', color='green')
    plt.yscale("log")
    # plt.ylim((0, 4000))
    # plt.xticks(x_pos, labels)
    plt.ylabel('Isomers')
    plt.xlabel('Family ID')
    plt.title('Isomeric Families per Set of Precursors.')

    iso_path = join(plot_path, 'Isomers_vs_Family_Id.png')
    plt.savefig(iso_path, dpi=600)
    plt.close()
    print('Save Isomers vs Family Id Plot at: "{}"'.format(iso_path))


def write_cytoscape_family_graph(family_dict):
    g_set    = set()
    pre_set  = set()
    g_lst = []
    core_set = set()
    all_set  = set()

    i = 0
    family_id = 0
    for prec_id, family_group in family_dict.items():
        family_pres = list(prec_id)
        
        for pre in prec_id:
            pre_set.add(pre)
        
        for f in range(len(family_pres) - 1, 0, -1):
            g_set.add((family_pres[f], family_pres[f-1], utils.string_diff(family_pres[f], family_pres[f-1]), family_id))
            
    #     for family in family_group:
    #         for line in family:
    #             core_set.add(line[1])
                
    #             prev = line[0]
    #             for j in range(len(line[2])):
    #                 frag = line[2][j]
    #                 all_set.add(frag)
    #                 g_set.add((prev, frag, string_diff(prev, frag), family_id))
    #                 prev = frag
                    
        family_id += 1

        g_lst = list(g_set)
    if g_lst:
        df = pd.DataFrame(g_lst, columns =['source', 'target', 'neutral_loss', 'family_id'])
        graph_path = join(file_path, "cytoscape_family_graph.csv")
        df.to_csv(graph_path, index=False)
        print('Wrote Cytoscape File at "{}"'.format(graph_path))


def family_parents_vs_oxygen_class(family_dict):
    prev = ""
    curr = ""
    parents = set()
    for family_id in family_dict:
        parents.add(family_id[-1])

    o_len = 20    
    o_count = [0] * o_len
    parents = set(parents)
    for pre in parents:
        if 'N' not in pre:
            o_count[utils.get_formula(pre)["O"]] += 1

    print("Distribution of Family Parents over Oxygen Class:\n{}".format(o_count))
    x_pos = range(3, 15)
    labels = range(3, 15)
    plt.bar(x_pos, o_count[3: 15], align='center', alpha=0.5)
    plt.xticks(x_pos, labels)
    plt.ylabel('Number of Family Parents')
    plt.xlabel('Oxygen Group')
    plt.title('Distribution of Family Parents over Oxygen Class.')
    family_path = join(plot_path, 'Family_Parents_vs_Oxygen_Class.png')
    plt.savefig(family_path, dpi=600)
    plt.close()
    print('Saved Family Parents vs Oxygen Class at: "{}"'.format(family_path))


def family_size_dist(family_dict):
    size = 0
    sizes = [0] * 10
    for family_id in family_dict:
        size = len(family_id)
        sizes[size] += 1

    print("Family size distribution:\n{}".format(sizes))
    x_pos = range(2, 7)
    y = sizes[2: 7]
    labels = range(2, 7)
    plt.bar(x_pos, y, align='center', alpha=0.5)
    plt.xticks(x_pos, labels)
    plt.ylabel('Number of Families')
    plt.xlabel('Family Size')
    plt.title('Distribution of Number of Families over Family Size')
    family_path = join(plot_path, 'Family_Size_Distribution.png')
    plt.savefig(family_path, dpi=600)
    plt.close()
    print('Saved Family Size Distribution at: "{}"'.format(family_path))


def family_dist_over_nl_seq(family_dict):
    fam_dict = defaultdict(int)
    for prec_id, family_group in family_dict.items():
        fam_set = set()
        for family in family_group:
            nl_str = ""
            first = True
            for row in family:
                if first:
                    first = False
                else:
                    nl_str += row[3][1:] + ", "
            fam_set.add(nl_str[:-2])
            
        for seq in fam_set:
            fam_dict[seq] += 1
    
    seqs = []
    seq_counts = []
    for seq, count in fam_dict.items():
        seqs.append(seq)
        seq_counts.append(count)

    seqs, seq_counts = list(zip(*sorted(zip(seqs, seq_counts), key=lambda x: x[1], reverse=True)))
    plt.bar(range(len(seqs)), seq_counts, align='center', color='green')
    plt.yscale("log")
    # plt.ylim((0, 4000))
    # plt.xticks(x_pos, labels)
    plt.ylabel('Number of Families')
    plt.xlabel('Neutral Loss Sequence ID')
    plt.title('Family Distribution over Neutral Loss Sequence')
    fam_path = join(plot_path, 'Neutral-Loss-Sequence-vs-Families.png')
    plt.savefig(fam_path, dpi=600)
    plt.close()
    print('Saved Family Count Distribution over Neutral Loss Sequence Plot at: "{}"'.format(fam_path))
    df = pd.DataFrame(list(zip(seqs, seq_counts)),
               columns =['Neutral Loss Sequence', 'Number of Families'])
    fam_path = join(file_path, 'Neutral-Loss-Sequence-vs-Families.csv')
    df.to_csv(fam_path)
    print('Saved Family Count Distribution over Neutral Loss Sequence csv at: "{}"'.format(fam_path))


def core_coverage(nominal_dict, pathway_dict, family_dict):
    total_cores = 0
    nom_core = {}
    for pre, cores in pathway_dict.items():
        nom_mz = nominal_dict[pre]
        for core in cores:        
            if nom_mz not in nom_core:
                nom_core[nom_mz] = set()
            nom_core[nom_mz].add(core)

    for _, core_set in nom_core.items():
        total_cores += len(core_set)

    nom_core = {}
    family_cores = 0
    for prec_id, family_group in family_dict.items():
        for family in family_group:
            for row in family:
                nom_mz = nominal_dict[row[0]]
                if nom_mz not in nom_core:
                    nom_core[nom_mz] = set()
                nom_core[nom_mz].add(row[1])
                
    for _, core_set in nom_core.items():
        family_cores += len(core_set)

    return total_cores, family_cores


def fragment_coverage(nominal_dict, pathway_dict, family_dict):
    pre_set = set()
    core_set = set()
    for pre, cores in pathway_dict.items():
        pre_set.add(pre)
        for core, frags in cores.items():
            core_set.add(core)
            
    total_frags = 0
    nom_frag = {}
    for pre, cores in pathway_dict.items():
        nom_mz = nominal_dict[pre]
        pre_frag_set = set()
        for core, paths in cores.items():
            for nl, path in paths.items():
                for frag in path["path"]:
    #                 if frag not in pre_set and frag not in core_set:
                    if nom_mz not in nom_frag:
                        nom_frag[nom_mz] = set()
                    nom_frag[nom_mz].add(frag)
                        
    for _, pre_frag_set in nom_frag.items():
        total_frags += len(pre_frag_set)

    family_frags = 0
    nom_frag = {}
    for prec_id, family_group in family_dict.items():
        pre_frag_set = set()
        for family in family_group:
            for row in family:
                nom_mz = nominal_dict[row[0]]
                for frag in row[2]:
    #                 if frag not in pre_set and frag not in core_set:
                    if nom_mz not in nom_frag:
                        nom_frag[nom_mz] = set()
                    nom_frag[nom_mz].add(frag)

    for _, pre_frag_set in nom_frag.items():
        family_frags += len(pre_frag_set)

    return total_frags, family_frags
    