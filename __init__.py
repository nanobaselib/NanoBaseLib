from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from nanobaselib.dataprep import dataprep
from nanobaselib.base_calling import base_calling
from nanobaselib.ploya_detection import ploya_detection
from nanobaselib.segment_align import segment_align
from nanobaselib.rna_mod_detection import rna_mod_detection

modules = [
    'dataprep', 'base_calling', 'ploya_detection', 'segment_align', 'rna_mod_detection'
]

__version__ = '0.1.0'


def main():
    parser = ArgumentParser(
        'nanobaselib',
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__)
    )

    subparsers = parser.add_subparsers(
        title='subcommands', description='valid commands',
        help='additional help', dest='command'
    )
    subparsers.required = True

    for module in modules:
        mod = globals()[module]
        p = subparsers.add_parser(module, parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    args = parser.parse_args()
    args.func(args)