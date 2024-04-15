#! /usr/bin/env python3

# AutoXTE
# Copyright 2023-2024 Nicholas Kuechel
# License Apache 2.0


import subprocess as sp
import os
import argparse
import pathlib as pl
from astroquery.heasarc import Heasarc
from astropy.table import Table
import glob
import shutil
import sys
import pkg_resources
import logging
import gzip

AUTOXTE = os.path.basename(sys.argv[0])
VERSION = pkg_resources.get_distribution('autoxte').version


def extract_gz(file) -> str:
    """
    Extracts .gz file
    """
    with gzip.open(file, "rb") as gz_in:
        fname = str(file).split(".gz")
        with open(fname[0], "wb") as orig_out:
            shutil.copyfileobj(gz_in, orig_out)
    os.remove(file)
    return f"{file} -> {fname[0]}"


def compress_gz(file) -> str:
    """
    Gzips a specified file and removes
    the original specified file after compression
    """
    with open(file, "rb") as f_in:
        with gzip.open(f"{file}.gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file)
    return f"{file} -> {file}.gz"


def query_obj(target=None):
    """
    Querys xtemaster @ heasarc for a specific target
    """
    if target is None:
        target = str(input("Target: "))
    heasarc = Heasarc()
    Heasarc.clear_cache()
    pca = heasarc.query_object(target, mission="xtemaster", resultmax=1000000)
    pca = Table(pca).to_pandas()
    pca["OBSID"] = pca["OBSID"].str.decode("utf-8")
    pca["PRNB"] = pca["PRNB"].str.decode("utf-8")
    pca = pca.loc[pca["EXPOSURE"] != 0]
    return pca


class Autoxte:
    def __init__(self, target=None):
        self.logger = logging.getLogger(AUTOXTE)
        self.logger.info("##############  Auto XTE  ##############\n")
        self.state = True
        self.table = query_obj(target)
        self.obs = []
        self.ras = []
        self.decs = []
        self.cNums = []
        self.prnbs = []
        self.base_dir = os.getcwd()

    def selection(self):
        """
        Interface for selecting the OBSIDs to be retrieved
        """
        cyc_aliass = ["cycle", "prnb"]
        while self.state is True:
            enter = str(input("autoXTE > ")).split(" ")
            if enter[0].lower() == "done":
                self.state = False
            elif enter[0].lower() == "sel":
                self.logger.info("Observations Selected:")
                for i in self.obs:
                    self.logger.info(i)
            elif enter[0].lower() == "back":
                del self.obs[-1]
                del self.ras[-1]
                del self.decs[-1]
                del self.cNums[-1]
                del self.prnbs[-1]
            elif enter[0].lower() in cyc_aliass:
                rows = self.table.loc[self.table["PRNB"] == str(enter[1])]
                self.logger.info(rows)
                for i in rows["OBSID"]:
                    self.obs.append(i)
                    self.logger.info(f"Added {i}")
                for i in rows["RA"]:
                    self.ras.append(i)
                for i in rows["DEC"]:
                    self.decs.append(i)
                for i in rows["PRNB"]:
                    self.prnbs.append(int(i))
                    prnb = str(i)
                    if int(prnb[0]) == 9:
                        cycles = [9, 10, 11, 12, 13, 14, 15, 16]
                        self.cNums.append(cycles[int(prnb[1])])
                    else:
                        self.cNums.append(int(prnb[0]))
            elif enter[0].lower() == "all":
                for i in self.table["OBSID"]:
                    self.obs.append(i)
                for i in self.table["RA"]:
                    self.ras.append(i)
                for i in self.table["DEC"]:
                    self.decs.append(i)
                for i in self.table["PRNB"]:
                    self.prnbs.append(i)
                    prnb = str(i)
                    if int(prnb[0]) == 9:
                        cycles = [9, 10, 11, 12, 13, 14, 15, 16]
                        self.cNums.append(cycles[int(prnb[1])])
                    else:
                        self.cNums.append(int(prnb[0]))
            elif enter[0].lower() == "exit":
                exit()
            else:
                try:
                    self.obs.append(enter[0])
                    row = self.table.loc[self.table["OBSID"] == enter[0]]
                    dt = row["TIME"]
                    row.reset_index(drop=True, inplace=True)
                    dt.reset_index(drop=True, inplace=True)
                    self.ras.append(row["RA"][0])
                    self.decs.append(row["DEC"][0])
                    self.prnbs.append(row["PRNB"][0])
                    prnb = str(row["PRNB"][0])
                    if int(prnb[0]) == 9:
                        cycles = [9, 10, 11, 12, 13, 14, 15, 16]
                        self.cNums.append(cycles[int(prnb[1])])
                    else:
                        self.cNums.append(int(prnb[0]))
                    self.ras.append(row["RA"][0])
                    self.decs.append(row["DEC"][0])
                except KeyError:
                    self.logger.info("Unknown Entry")

    def pull(self, index):
        """
        Retrieves an XTE dataset for a specified obsid
        """
        obsid = self.obs[index]
        downurl = "https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive/AO"
        self.logger.info("---------------------------------------------------")
        self.logger.info(f"          Prosessing OBSID: {obsid}")
        self.logger.info("---------------------------------------------------")
        self.logger.info("Downloading 00 data...")
        fullurl = (
            f"{downurl}{self.cNums[index]}//" +
            f"P{self.prnbs[index]}/{obsid}/"
        )
        endargs = "--show-progress --progress=bar:force"
        downcommand = (
            "wget -q -nH --no-check-certificate --cut-dirs=5 "
            + "-r -l0 -c -N -np -R 'index*'"
            + f" -erobots=off --retr-symlinks {fullurl} {endargs}"
        )
        return sp.call(downcommand, shell=True)

    def barycenter(self, index):
        """
        Calls the HEASoft's barycorr tool to perform
        a barycenter correction to .evt data
        """
        obsid = self.obs[index]
        obsbasepath = f"P{self.prnbs[index]}/{self.obs[index]}"
        try:
            orbfile = glob.glob(f"{self.base_dir}/{obsbasepath}/orbit/FP*")[0]
        except IndexError:
            self.logger.info("Orbitfile not found")
            return False
        self.logger.info(f"Orbitfile: {orbfile}")
        pcapath = pl.Path(f"{obsbasepath}/pca/").resolve()
        os.chdir(pcapath)
        gzevts = glob.glob("*.evt.gz")
        for i in gzevts:
            extract_gz(i)
        evt_files = glob.glob("*evt")
        for n, i in enumerate(evt_files, start=1):
            evttype = i.split("_")[0]
            sp.call(
                f"barycorr infile={i} outfile=bc{evttype}_{obsid}_{n}.evt"
                + f" orbitfiles={orbfile} refframe=ICRS"
                + f" ra={self.ras[index]} dec={self.decs[index]}",
                shell=True,
            )
            comp_org = compress_gz(i)
            self.logger.info(comp_org)
        os.chdir(self.base_dir)
        return True

    def remove_bloatdirs(self, index) -> None:
        """
        Removed unnecessary dirs in pulled XTE datasets
        as defined at
        https://heasarc.gsfc.nasa.gov/docs/xte/start_guide.html,
        section 4.5
        """
        obsbasepath = f"P{self.prnbs[index]}/{self.obs[index]}"
        bloatdirs = ["ace", "acs", "clock", "eds", "fds", "gsace",
                     "ifog", "ipsdu", "pse", "spsdu"]
        for i in bloatdirs:
            shutil.rmtree(f"{obsbasepath}/{i}/")

    def pull_reduce(self, bc) -> None:
        """
        Pulls all selected OBSIDs
        """
        for i in self.obs:
            index = self.obs.index(i)
            self.pull(index)
            if bc is True:
                self.barycenter(index)
            self.remove_bloatdirs(index)


def cli(args=None):
    """
    Command line interface argument parsing
    """
    p = argparse.ArgumentParser(
            description="A program for piplining XTE data" +
                        "retrieval and barycenter correcting"
    )

    p.add_argument(
        "-src",
        "--source",
        help="Set the src from the command line",
        type=str,
        default=None,
    )

    p.add_argument(
        "-bc",
        "--bc",
        help="Engages option for a barycenter correction" +
             "in data reduction procedure",
        action="store_true",
        default=None,
    )
    argp = p.parse_args(args)
    return argp


def main():
    level = logging.INFO
    logging.basicConfig(stream=sys.stdout,
                        level=level,
                        format="")
    args = cli()
    xte = Autoxte(target=args.source)
    xte.selection()
    xte.pull_reduce(args.bc)
