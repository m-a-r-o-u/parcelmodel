import sys
import yaml
from setup import model_init
from logger import logger_factory

def logger_config(config):
  return config['output']

def main():
    config = yaml.load(sys.stdin)
    model = model_init(config['initial_conditions'])
    with logger_factory(logger_config(config)) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
