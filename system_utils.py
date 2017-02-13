import os

def check_make_directory(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)
