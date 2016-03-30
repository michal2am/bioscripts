# python 3
# functions for file handling
# michaladammichalowski@gmail.com
# 30.03.15 - creation
# EXAMPLE CALL:

import os


def find_files(directory, extension):
    """
    :param directory: starting directory
    :param extension: files to look for
    :return: list of relative paths to files
    """
    file_list = []
    for dir_, _, files in os.walk(directory):
        for file in files:
            if file.endswith(extension):
                relDir = os.path.relpath(dir_, directory)
                file_list.append(os.path.join(relDir, file))
    return file_list

def create_outnames(files, postfix):
    """
    :param files: list of files to change extension
    :param postfixes: new extension
    :return:list of relative paths to files with changes extension
    """
    return [os.path.splitext(fil)[0] + postfix for fil in files]
