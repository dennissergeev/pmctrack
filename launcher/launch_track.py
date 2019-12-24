#!/usr/bin/env python
"""
Submit tracking jobs to a PBS cluster queue
"""
# from datetime import datetime
import json
import os
import subprocess as sb

aux_dict = dict(dt_start="2010_10_01_0000", dt_end="2011_04_30_2300")
winter = aux_dict["dt_start"][:4] + "_" + aux_dict["dt_end"][:4]

DATADIR = os.path.join(os.getenv("PROJ",), "reanalysis", "{dataset}")
OUTDIR = os.path.join(
    os.getenv("PROJ"),
    "pmctrack",
    "output",
    "{dataset}",
    "{datestamp}run{run:03d}",
    "{winter}",
)
STORE_DIR = "configs"
datasets = ["era5", "interim"]
walltimes = ["4:0:0", "0:15:0"]
DEF_CONF_FNAME = os.path.join(STORE_DIR, "dummy_settings.conf")
CONF_MASK = os.path.join(
    STORE_DIR, ("{dataset}_{datestamp}run{run:03d}_{winter}" "_settings.conf")
)

LOG_DIR = os.path.join(os.getenv("PROJ"), "pmctrack", "logs")
EXE = os.path.join(os.getenv("HOME"), "pmctrack", "track.x")
# cd $PBS_O_WORKDIR
CMD_MASK = (
    "qsub"
    " -N pmctrack_{dataset}_{datestamp}run{run:03d}_{winter}"
    " -q shared"
    " -l select=1"
    " -l walltime={walltime}"
    " -o {log_dir}"
    " -e {log_dir}"
    " -- {exe} {conf}"
)

# today = datetime.now()
# datestamp = '{dt:%Y%m%d}_'.format(dt=today)
datestamp = ""

with open("runs_grid.json", "r") as jf:
    runs_list = json.load(jf)


for dataset, walltime in zip(datasets, walltimes):
    with open(DEF_CONF_FNAME, "r") as fp:
        def_conf = fp.readlines()

    aux_dict["datadir"] = '"{}"'.format(DATADIR.format(dataset=dataset))
    aux_dict["fname_lvl"] = '"{dataset}.an.pl.%YYYY.%MM.%VAR.nc"'.format(
        dataset=dataset
    )
    aux_dict["fname_sfc"] = '"{dataset}.an.sfc.%YYYY.%MM.%VAR.nc"'.format(
        dataset=dataset
    )

    for run, params_dict in enumerate(runs_list):
        aux_dict["outdir"] = '"{}"'.format(
            OUTDIR.format(dataset=dataset, run=run, datestamp=datestamp, winter=winter)
        )
        target_conf = []
        for line in def_conf:
            new_l = line.rstrip("\n")
            if not new_l.startswith("#"):
                old_param = new_l.split("=")[0]
                d = params_dict.copy()
                d.update(aux_dict)
                for param, new_val in d.items():
                    if old_param == param:
                        new_l = "{}={}".format(param, new_val)
            target_conf.append(new_l)

        conf_file = os.path.abspath(
            CONF_MASK.format(
                dataset=dataset, datestamp=datestamp, winter=winter, run=run
            )
        )
        with open(conf_file, "w") as fp:
            for line in target_conf:
                fp.write(line + "\n")

        subp_args = CMD_MASK.format(
            dataset=dataset,
            run=run,
            winter=winter,
            datestamp=datestamp,
            walltime=walltime,
            log_dir=LOG_DIR,
            exe=EXE,
            conf=conf_file,
        ).split()
        # ['cp {} {}'.format(conf_file, aux_dict['outdir'])])
        print(" ".join(subp_args))
        sb.call(subp_args)
