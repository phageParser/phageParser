import pickle


def find_incomplete_records(gen_dict):
    incomp_list = []
    for loc_id in gen_dict:
        locus = gen_dict[loc_id]
        if set(locus.keys()) != {'Spacers', 'Start', 'RepeatSeq', 'Stop'}:
            incomp_list.append(loc_id)
    return incomp_list


def find_size_offsets(gen_dict):
    offset_list = []
    for loc_id in gen_dict:
        locus = gen_dict[loc_id]
        pos_size = int(locus['Stop']) - int(locus['Start'])
        spacer_size = sum(len(s) for s in locus['Spacers'].values())
        repeat_size = len(locus['RepeatSeq']) * (len(locus['Spacers']) + 1)
        size_diff = pos_size - spacer_size - repeat_size
        offset_list.append((size_diff, loc_id))
    return offset_list


def del_keys(dictionary, key_list):
    for key in key_list:
        if key in dictionary:
            del dictionary[key]
    return dictionary


def prune_dict(gen_dict):
    # remove loci with incomplete fields
    incomp_list = find_incomplete_records(gen_dict)
    gen_dict = del_keys(gen_dict, incomp_list)
    # remove loci with offsets greater than 0
    offset_list = find_size_offsets(gen_dict)
    locitodel = [x[1] for x in offset_list if x[0] > 0]
    gen_dict = del_keys(gen_dict, locitodel)
    return gen_dict


if __name__ == '__main__':
    with open('gendict.pickle', 'rb') as f:
        gen_dict = pickle.load(f)
    gen_dict = prune_dict(gen_dict)
    with open('gendictpruned.pickle', 'wb') as f:
        pickle.dump(gen_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
