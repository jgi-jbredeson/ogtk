
def count_items(items):
    count = {}
    for item in items:
        if item in count:
            count[item] += 1
        else:
            count[item] = 1
    return count


def index_list(list_, start=0):
    return dict(zip(list_, range(start,start+len(list_))))

