from src import config


def overlap_size(path1, path2):
    len1 = len(path1)
    len2 = len(path2)
    
    min_len = min(len1, len2)
    path1 = path1[:min_len]
    path2 = path2[:min_len]
    
    overlap_count = 0
    for frag1, frag2 in zip(path1, path2):
        if frag1 == frag2:
            overlap_count += 1
        else:
            break
    
    return overlap_count


def get_key(row, overlap_len, key_offset):
    precursor = row["Precursor"]
    
    if len(row["Pathway"]) >= overlap_len:
        path = row["Pathway"][:overlap_len]
        short_key = False
    else:
        path = row["Pathway"]
        short_key = True
    
    if key_offset < 0:
        path = path[:key_offset]
    
    key = (precursor, tuple(path))
    return key, short_key


def in_graph(key, graph):
    if key in graph:
        return True, key, 0
    for i in range(len(key[1]) - 1):
        short_key = (key[0], key[1][:-(i+1)])
        if (short_key in graph) and (graph[short_key]["short-key"] == True):
            return True, short_key, -(i+1)
    return False, key, 0


def add_node(row, graph, overlap_len, key_offset, roots, is_root):
    key, short_key = get_key(row, overlap_len, key_offset)
    
    if key not in graph:
        graph[key] = {"core":row["Core-Fragment"], "pathway":row["Pathway"], "fnl":row["ID"].split(' ')[0], "short-key":short_key, "edges":set()}
        if is_root:
            roots.add(key)
    return key


def add_edge(key, parent_key, original_pathway, graph):
    assert parent_key in graph
    graph[parent_key]["edges"].add((key, original_pathway))


def get_path_forest(pathway_dict_df):
    path_forest = {}
    roots = set()

    fnl = ["1CO2", "1CH4O", "1H2O", "1CO"]

    for index, row in pathway_dict_df.iterrows():
        neutral_diff = row["ID"].split(' ')[0]
        if neutral_diff not in fnl:
            add_node(row, path_forest, config.overlap_len, key_offset=0, roots=roots, is_root=True)
            continue
            
        parent_prec = row["Pathway"][0]
        if len(row["Pathway"][1:]) >= config.overlap_len:
            parent_path = row["Pathway"][1:config.overlap_len+1]
        else:
            parent_path = row["Pathway"][1:]
            
        parent_key = (parent_prec, tuple(parent_path))
        parent_exists, parent_key, key_offset = in_graph(parent_key, path_forest)
        if parent_exists:
            key = add_node(row, path_forest, config.overlap_len, key_offset, roots, is_root=False)
            add_edge(key, parent_key, tuple(row["Pathway"]), path_forest)
        else:
            add_node(row, path_forest, config.overlap_len, key_offset=0, roots=roots, is_root=True)

    return roots, path_forest


def populate_family_dict(path, path_forest, family_dict):
    family = []
    prec_list = []
    for node in path:
        row = path_forest[node[0]]
        out_row = [node[0][0], row["core"], list(node[1]), row["fnl"], node[0][1], row["short-key"]]
        prec_list.append(node[0][0])
        family.append(out_row)
    prec_id = tuple(prec_list)
    if prec_id in family_dict:
        family_dict[prec_id].append(family)
    else:
        family_dict[prec_id] = [family]
    # if prec_id in family_dict:
    #     family_dict[prec_id] += 1
    # else:
    #     family_dict[prec_id] = 1


def combine_families(roots, path_forest):
    # this code is combining families with same set of precurosrs
    family_counter = 0
    family_dict = {}
    for root in roots:
        if not path_forest[root]["edges"]:
            continue
        stack = [((root, tuple(path_forest[root]["pathway"])), None)]
        path = []
        while stack:
            node, _ = stack.pop()
            path.append(node)
            children = path_forest[node[0]]["edges"]
            if children:
                for child in children:
                    stack.append((child, node))
            else:
                family_counter += 1
                populate_family_dict(path, path_forest, family_dict)
                path.pop()
                while stack and path[-1] != stack[-1][1]:
                    path.pop()
    
    return family_dict