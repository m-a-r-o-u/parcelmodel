def modify_config_tree(inp, year):
    if isinstance(inp, list):
        return [modify_config_tree(i, year) for i in inp]
    if isinstance(inp, dict):
        return {k:modify_config_tree(v, year) for k,v in inp.iteritems()}
    if isinstance(inp, str):
        return inp.replace('time_stamp', year)
    return inp

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
    return config

def logger_config(config):
    modified_config = set_timestamp(config)
    modified_config = set_global_path(modified_config)
    return modified_config['output']['logger']

