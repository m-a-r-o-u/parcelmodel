import sys
import yaml
from setup import model_init
from logger import logger_factory

def return_timestamp():
  '''returns the current time stamp'''
  import time, datetime

  tt=time.time()
  dt=datetime.datetime.fromtimestamp(tt).strftime('%Y-%m-%dT%H:%M:%S')
  return dt

def set_timestamp(config):
    path_key = 'file_path'
    current_time_stamp = return_timestamp()
    if path_key in config['output']:
        config['output'][path_key] = config['output'][path_key].replace('time_stamp', current_time_stamp)
    return config

def set_global_path(config):
    path_key = 'file_path'
    output_path = None
    if path_key in config['output']:
        output_path =  config['output'][path_key]
    if not output_path == None:
        for c in config['output']['logger']:
          if path_key in c:
            c[path_key] = output_path
    return config

def logger_config(config):
    config = set_timestamp(config)
    config = set_global_path(config)
    return config['output']['logger']

def main():
    config = yaml.load(sys.stdin)
    model = model_init(config['initial_conditions'])
    with logger_factory(logger_config(config)) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
