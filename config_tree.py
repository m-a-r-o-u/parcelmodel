def modify_config_tree(inp, match, f, path=None):
    if path is None:
        path = []
    if match(path, inp):
        return f(inp)
    if isinstance(inp, list):
        return [modify_config_tree(element, match, f, path+[i]) for i, element in enumerate(inp)]
    if isinstance(inp, dict):
        return {k:modify_config_tree(v, match, f, path+[k]) for k,v in inp.iteritems()}
    return inp

def replace_string_in_tree(inp, a, b):
    return modify_config_tree(inp, lambda _path, element: isinstance(element, str), lambda x: x.replace(a, b))

def set_dict_item_in_tree(inp, key, value, path_prefix=None):
    if path_prefix is None:
        path_prefix = []
    def match(path, element):
        return isinstance(element, dict) and path[:len(path_prefix)]==path_prefix
    def modifier(x):
        x = x.copy()
        x[key] = value
        return x
    return modify_config_tree(inp, match, modifier)

def modify_item_in_config_tree(inp, key, value):
    if isinstance(inp, list):
        return [modify_item_in_config_tree(i, key, value) for i in inp]
    return {k:modify_one_item(k, v, key, value) for k,v in inp.iteritems()}

def modify_one_item(k, v, k_ref, v_ref):
    if k == k_ref:
        return v_ref
    else:
        return v

def return_timestamp():
  '''returns the current time stamp'''
  import time, datetime
  tt=time.time()
  dt=datetime.datetime.fromtimestamp(tt).strftime('%Y-%m-%dT%H:%M:%S')
  return dt

def set_timestamp(config_tree):
    current_timestamp = return_timestamp()
    return replace_string_in_tree(config_tree, 'time_stamp', current_timestamp)

def set_global_path(config):
    key = 'file_path'
    if key in config['globals']:
        value = config['globals'][key]
        config = set_dict_item_in_tree(config, key, value, ['output', 'logger'])
    return config

def logger_config(config):
    modified_config = set_timestamp(config)
    modified_config = set_global_path(modified_config)
    return modified_config['output']['logger']
