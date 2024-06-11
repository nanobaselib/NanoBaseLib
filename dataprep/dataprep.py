import os
import re
from shutil import rmtree
from zipfile import ZipFile
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from bonito.util import __data__, __models__
from bonito.cli.convert import main as convert
from bonito.cli.convert import argparser as cargparser

import requests
from tqdm import tqdm


def main(args):
    """
    Download models and training sets
    """
    if args.models or args.all:

        if args.show:
            print("[available models]")
            for model in models:
                print(f" - {model}")
        else:
            print("[downloading models]")
            for model in models.values():
                File(__models__, model, args.force).download()
    if args.training or args.all:
        print("[downloading training data]")
        for train in training:
            File(__data__, train, args.force).download()


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--all', action='store_true')
    group.add_argument('--models', action='store_true')
    group.add_argument('--training', action='store_true')
    parser.add_argument('--list', '--show', dest='show', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    return parser