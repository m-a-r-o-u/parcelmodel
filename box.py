import numpy as np
import boxmodel_functions as bf
import boxmodel_constants as bc
from model import Model
from logger import Logger, FinalStateLogger, PlotTLogger, PlotQVLogger, MultiLogger

def logger_init(final_only=True):
    if final_only==0:
        return FinalStateLogger()
    elif final_only==1:
        return Logger()
    elif final_only==2:
        return PlotTLogger()
    elif final_only==3:
        return PlotQVLogger()
    elif final_only==4:
        return MultiLogger([PlotTLogger(), PlotQVLogger()])

def model_init():
    return Model()

def main():
    logger = logger_init(2)
    model = model_init()
    model.run(logger)
    logger.finalize()

if __name__ == "__main__":
    main()
