import argparse


def argparse_init():
    parser = argparse.ArgumentParser(description='process yaml-config')
    parser.add_argument('yaml', help='yaml file needed for configuration')
    args = parser.parse_args()
    return args
