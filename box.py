import sys
import yaml
from model import Model
from logger import logger_factory

def model_init():
    return Model()

def main():
    config = yaml.load(sys.stdin)
    model = model_init()
    with logger_factory(config['output']) as logger:
        model.run(logger)

if __name__ == "__main__":
    main()
