import yaml
from .setup import model_init
from .logger import logger_factory
from .config_tree import logger_config
from .my_argparse import argparse_init
from . import executer


def main():
    args = argparse_init()
    config = yaml.load(open(args.yaml), Loader=yaml.FullLoader)
    model = model_init(config['initial_conditions'],
                       executer.DISPATCHER[config['globals']['executer']])
    with logger_factory(logger_config(config)) as logger:
        logger.inform(config['initial_conditions'])
        model.run(logger)


if __name__ == "__main__":
    main()
