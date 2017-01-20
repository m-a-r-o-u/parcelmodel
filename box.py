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

def logger_config(config):
  return config['output']['logger']

def main():
    config = yaml.load(sys.stdin)
    model = model_init(config['initial_conditions'])
    with logger_factory(logger_config(config)) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
