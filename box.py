import sys
import yaml
from setup import model_init
from logger import logger_factory
from config_tree import logger_config
from my_argparse import argparse_init

def main():
    args = argparse_init()
    config = yaml.load( open(args.yaml) )
    model = model_init(config['initial_conditions'])
    with logger_factory(logger_config(config)) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
