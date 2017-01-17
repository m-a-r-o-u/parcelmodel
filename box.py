import sys
import yaml
from setup import model_init
from logger import logger_factory

def main():
    config = yaml.load(sys.stdin)
    model = model_init(config['initial_conditions'])
    with logger_factory(config['output']) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
