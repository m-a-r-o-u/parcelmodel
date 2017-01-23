def modify_config_tree(inp, year):
    if isinstance(inp, list):
        return [modify_config_tree(i, year) for i in inp]
    if isinstance(inp, dict):
        return {k:modify_config_tree(v, year) for k,v in inp.iteritems()}
    if isinstance(inp, str):
        return inp.replace('time_stamp', year)
    return inp

def modify_item_in_config_tree(inp, key, value):
    if isinstance(inp, list):
        return [modify_item_in_config_tree(i, key, value) for i in inp]
    mod_inp = inp
    mod_inp[key] = value
    return mod_inp

def return_timestamp():
  '''returns the current time stamp'''
  import time, datetime
  tt=time.time()
  dt=datetime.datetime.fromtimestamp(tt).strftime('%Y-%m-%dT%H:%M:%S')
  return dt

def set_timestamp(config):
    current_timestamp = return_timestamp()
    return modify_config_tree(config, current_timestamp)

def set_global_path(config):
    key = 'file_path'
    if key in config['globals']:
        value = config['globals']['file_path']
        modify_item_in_config_tree(config['output']['logger'], key, value)
    return config

def logger_config(config):
    modified_config = set_timestamp(config)
    modified_config = set_global_path(modified_config)
    return modified_config['output']['logger']

