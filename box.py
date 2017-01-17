import sys
import yaml
from model import Model
from logger import logger_factory

def model_init():
    return Model()

def main():
    config = yaml.load(sys.stdin)
    logger = logger_factory(config['output'])
    model = model_init()
    model.run(logger)
    logger.finalize()

if __name__ == "__main__":
    main()
