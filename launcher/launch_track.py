#!/usr/bin/env python
"""
Submit tracking jobs to a PBS cluster queue
"""
import json
import os
import subprocess as sb

aux_dict = dict(
    dt_start='2010_10_01_0000',
    dt_end='2011_04_30_2300'
)
winter = '2010_2011'

DATADIR = os.path.join(os.getenv('PROJ'), 'reanalysis', '{dataset}')
OUTDIR = os.path.join(os.getenv('PROJ'), 'pmctrack', 'output', '{dataset}',
                      '{run}', '{winter}')
STORE_DIR = 'configs'
datasets = ['era5', 'interim']
DEF_CONF_FNAME = os.path.join(STORE_DIR, 'settings.conf')
CONF_MASK = os.path.join(STORE_DIR, '{dataset}_{run}_settings.conf')


with open('runs_grid.json', 'r') as jf:
    runs_dict = json.load(jf)


for dataset in datasets:
    with open(DEF_CONF_FNAME, 'r') as fp:
        def_conf = fp.readlines()

    aux_dict['datadir'] = DATADIR.format(dataset=dataset)
    aux_dict['prefix_lvl'] = '{dataset}.an.pl.'.format(dataset=dataset)
    aux_dict['prefix_sfc'] = '{dataset}.an.sfc.'.format(dataset=dataset)

    for run, params_dict in runs_dict.items():
        aux_dict['outdir'] = OUTDIR.format(dataset=dataset, run=run,
                                           winter=winter)
        target_conf = []
        for line in def_conf:
            new_l = line.rstrip('\n')
            if not line.startswith('#'):
                old_param, = line.split('=')
                d = dict(**params_dict, **aux_dict)
                for param, new_val in d.items():
                    if old_param == param:
                        new_l = '{}={}'.format(param, new_val)
            target_conf.append(new_l)

        conf_file = CONF_MASK.format(dataset=dataset, run=run)
        with open(conf_file, 'w') as fp:
            for line in target_conf:
                fp.write(line+'\n')

        subp_args = ([])
        print(subp_args)
