# -*- coding: utf-8 -*-
"""
Miscellaneous utilities
"""
import numpy as np
from path import Path
import subprocess
#
# TODO: doc strings!
#


class PMCSettings:
    def __init__(self, fname_path):
        with fname_path.open('r') as f:
            self.conf = [line for line in f.read().split('\n')
                         if not line.startswith('#') and len(line) > 0]
        for line in self.conf:
            try:
                exec(line, None, self.__dict__)
            except SyntaxError:
                k, v = line.split('=')
                self.__dict__.update({k: str(v)})

    def __repr__(self):
        return ('Settings used for '
                f'PMC tracking algorithm ({len(self.conf)})')

    def __str__(self):
        summary = '\n'.join(self.conf)
        return f'Settings used for PMC tracking algorithm:\n\n{summary}'


def unit_format(value, unit='1'):
    if value == 1:
        if unit == '1':
            string = ''
        else:
            string = f'${str(unit).replace(" ", "$ $")}$'
            if '%' in string:
                string = string.replace("%", "\%")
    else:
        exp = np.floor(np.log10(value))
        base = value / 10**exp
        if exp == 0 or exp == 1:
            string = r'${0}$'.format(value)
        elif exp == -1:
            string = r'${0:0.1f}$'.format(value)
        else:
            if int(base) == 1:
                string = r'$10^{{{0:d}}}$'.format(int(exp))
            else:
                string = r'${0:d}\times10^{{{1:d}}}$'.format(int(base),
                                                             int(exp))
        if not unit == '1':
            string += r' ${}$'.format(unit)
    return string


def merge_to_gif(input_images, gifname, delay=20, resize=100, output_dir=None):
    if isinstance(input_images, str):
        input_path = Path(input_images)
    if output_dir is None:
        output_dir = input_path.dirname()
    else:
        output_dir = Path(output_dir)
    subp_args = (['convert', '-resize', str(resize)+'%',
                             '-delay', str(delay),
                             '-dispose', 'Background',
                             '+page', input_images,
                             '-loop', '0', output_dir / gifname])
    print('Generating .gif: {0}'.format(' '.join(subp_args)))
    subprocess.check_call(subp_args)


def _make_dir(d, force=False):
    if d.exists():
        if force:
            d.rmtree_p()
            d.makedirs()
        # else:
        #     warnings.warn('Directory {d} exists!'.format(d=d))
    else:
        d.makedirs()
