def group_by(iterable, key=lambda x: (x, x)):

    items = [key(el) for el in iterable]
    groups = {k: [] for k, v in items}
    [groups[k].append(v) for k, v in items]
    return groups


def remove_items(dictionary, cond):

    removed_items = []
    for item in list(dictionary.items()):
        if cond(item):
            dictionary.pop(item[0])
            removed_items.append(item[1])
    return removed_items