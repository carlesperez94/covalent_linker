import os

def create_file(content, filename):
    with open(filename, "w") as out:
        out.write(content)

def get_directory_from_filepath(path):
    return os.path.dirname(path)
