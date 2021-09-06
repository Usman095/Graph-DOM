import re
from collections import OrderedDict

import config

# get mass of a molecule given in dictionary form.
def get_mass(dic_formula):
    mass = 0
    for item in dic_formula.items():
        mass += config.elements[item[0]] * item[1]
    return mass


# get formula dictionary from string
def get_formula(formula, mol=config.atom_list):
    
    res = OrderedDict.fromkeys(mol, 0)  
    formula_list = re.findall(r'([A-Z][a-z]*)(-?\d*)', formula)
    for key, value in formula_list:
        res[key] = int(value) if value else 1
    
    return res


# calculate the difference between two formulas provided as dictionary with keys C, H, O
# formula_a - formula_b
def formula_diff(formula_a, formula_b, mol=config.atom_list):
    temp_a = OrderedDict()
    temp_b = OrderedDict()
    for atom in mol:
        temp_a[atom] = formula_a[atom] if atom in formula_a else 0
        temp_b[atom] = formula_b[atom] if atom in formula_b else 0
    
    diff_mol = OrderedDict()
    for atom in mol:
        diff_mol[atom] = temp_a[atom] - temp_b[atom]
    
    return diff_mol


# FIXME: standardize the formula. Add mol parameter.
def get_string_from_formula(formula):
    str_lst = []
    for element, count in formula.items():
        if count > 0:
            str_lst.append(str(element))
            str_lst.append(str(count))
    return ''.join([elem for elem in str_lst])


def string_diff(a, b, mol=config.atom_list):
    a_f = get_formula(a, mol)
    b_f = get_formula(b, mol)
    diff_f = formula_diff(a_f, b_f, mol)
    return get_string_from_formula(diff_f)


#get the atom count in mol (single letter atoms only C, H, O etc.)
def get_count(atom, mol):
    
    if atom not in mol:
        return 0
    
    count = mol.split(atom)[1]
    if len(count) == 0 or count[0].isalpha(): #if it's last element or the next element is a character
        return 1
    
    result = re.findall('\d+|$', mol.split(atom)[1])[0]
    return int(result) if result else 0


#get the core fragment int from a precursor formula dict and a dict of neutral loss counts
def get_core_int(nl_dict, precursor_formula):
    core_frag = []
    for element, count in precursor_formula.items():
        element_tot_count = 0
        for nl, _ in config.neutral_losses.items():
            element_nl_count = get_count(element, nl)
            element_tot_count += element_nl_count * nl_dict[nl]
            #if index == 0:
            #    print('Number of {}\'s in neutral loss {} = {} * {}'.format(element, nl, element_nl_count, row[nl]))
        
        # FIXME: "count - element_tot_count > 0" and "element_tot_count < count" is the same.
        if count - element_tot_count > 0 and element_tot_count < count:
            cho = str(count - element_tot_count)
            core_frag.append(cho if len(cho) == 2 else '0{}'.format(cho))
    
    return int(''.join(core_frag))


#get the core fragment string from a precursor formula dict and a dict of neutral loss counts
def get_core_string(nl_dict, precursor_formula):
    core_frag = []
    at_least_one = False
    for element, count in precursor_formula.items():
        element_tot_count = 0
        all_losses = OrderedDict()
        all_losses.update(config.neutral_losses)
        all_losses.update(config.alt_losses)
        for nl, _ in all_losses.items():
            element_nl_count = get_count(element, nl)
            element_tot_count += element_nl_count * nl_dict[nl]
            #if index == 0:
            #    print('Number of {}\'s in neutral loss {} = {} * {}'.format(element, nl, element_nl_count, row[nl]))
        
        
        # FIXME: "count - element_tot_count > 0" and "element_tot_count < count" is the same.
        if element_tot_count < count:
            at_least_one = True
            cho = str(count - element_tot_count)
            core_frag.append(element)
            core_frag.append(cho)
        
        if element_tot_count > count:
            raise ValueError(f"For precursor {get_string_from_formula(precursor_formula)} and pathway {nl_dict}, \
                             total count of element {element} is too large: {element_tot_count}.")
    if not at_least_one:
        raise ValueError(f"core_fragment return null")
    
    return ''.join(core_frag)


"""get the core fragment string from a precursor formula dict and a list of tuples of neutral loss counts"""
def get_core_string_2(nl_list, precursor_formula):
    """convert the list of tuples to dictionary and then call the get_core_string function."""
    keys = []
    keys.extend(config.neutral_losses.keys())
    keys.extend(config.alt_losses.keys())
    nl_dict = OrderedDict.fromkeys(keys, 0)
    for nl, count in nl_list:
        if nl in nl_dict:
            nl_dict[nl] += count
        else:
            raise KeyError(f"{nl} is not a valid neutral loss. \
                           Only the following neutral losses are accepted at the moment: {config.neutral_losses.keys()}")
            
    return get_core_string(nl_dict, precursor_formula)


#get a sequence of intermediate fragments from the sequence of neutral losses and precursor ion (dictionary).
def get_fragment_seq_from_nloss(nloss_seq, precursor_formula):
    current_fragment_formula = precursor_formula.copy()
    fragment_seq = []
    
    for nloss, count in nloss_seq:
        assert count > 0
        for _ in range(count):
            nloss_formula = get_formula(nloss)
            current_fragment_formula = formula_diff(current_fragment_formula, nloss_formula)
            
        current_fragment = get_string_from_formula(current_fragment_formula)
        fragment_seq.append(current_fragment)
    
    return fragment_seq


def is_path_valid(path, precursor_formula): 
    keys = []
    keys.extend(config.neutral_losses.keys())
    keys.extend(config.alt_losses.keys())
    nl_dict = OrderedDict.fromkeys(keys, 0)
    for nl, count in path:
        if nl in nl_dict:
            nl_dict[nl] += count
        else:
            raise KeyError(f"{nl} is not a valid neutral loss. \
                           Only the following neutral losses are accepted at the moment: {nl_dict.keys()}")
    is_valid = True
    for element, count in precursor_formula.items():
        element_tot_count = 0
        all_losses = OrderedDict()
        all_losses.update(config.neutral_losses)
        all_losses.update(config.alt_losses)
        for nl, _ in all_losses.items():
            element_nl_count = get_count(element, nl)
            element_tot_count += element_nl_count * nl_dict[nl]
        if element_tot_count > count:
            is_valid = False
    return is_valid
