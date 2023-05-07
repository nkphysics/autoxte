#! /usr/bin/env python3

# AutoXTE
# Copyright 2023 Nicholas Kuechel
# License Apache 2.0


import subprocess as sp
import os
import argparse
import pandas
import pathlib as pl
from astroquery.heasarc import Heasarc
from astropy.table import Table
import glob
import shutil
import gzip


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
        print("##############  Auto XTE  ##############\n")
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
        while self.state is True:
            enter = str(input("autoXTE > "))
            if enter.lower() == "done":
                self.state = False
            elif enter == "version" or enter == "-v":
                print("autoXTE V4.1 Date Modified September 14, 2021")
            elif enter.lower() == "sel":
                print("Observations Selected:")
                for i in self.obs:
                    print(i)
            elif enter.lower() == "back":
                del self.obs[-1]
                del self.ras[-1]
                del self.decs[-1]
                del self.cNums[-1]
                del self.prnbs[-1]
            else:
                if len(str(enter)) == 5:
                    rows = self.table.loc[self.table["PRNB"] == int(enter)]
                    print(rows)
                    for i in rows["OBSID"]:
                        self.obs.append(i)
                    for i in rows["RA"]:
                        self.ras.append(i)
                    for i in rows["DEC"]:
                        self.decs.append(i)
                    for i in rows["PRNB"]:
                        self.prnbs.append(i)
                        prnb = str(i)
                        if int(prnb[0]) == 9:
                            cycles = [9, 10, 11, 12, 13, 14, 15, 16]
                            self.cNums.append(cycles[int(prnb[1])])
                        else:
                            self.cNums.append(int(prnb[0]))
                elif enter == "all":
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
                else:
                    self.obs.append(enter)
                    row = self.table.loc[self.table["OBSID"] == enter]
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

    def pull(self, index):
        """
        Retrieves an XTE dataset for a specified obsid
        """
        obsid = self.obs[index]
        downurl = "https://heasarc.gsfc.nasa.gov/FTP/xte/data/archive/AO"
        print("-----------------------------------------------------")
        print(f"           Prosessing OBSID: {obsid}")
        print("-----------------------------------------------------")
        print("Downloading 00 data...")
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
        orbfile = glob.glob(f"{self.base_dir}/{obsbasepath}/orbit/FP*")[0]
        print(f"Orbitfile: {orbfile}")
        pcapath = pl.Path(f"{obsbasepath}/pca/").resolve()
        os.chdir(pcapath)
        gzevts = glob.glob("*.evt.gz")
        for i in gzevts:
            extract_gz(i)
        evt_files = glob.glob("*evt")
        ccn = 1
        for i in evt_files:
            evttype = i.split("_")[0]
            sp.call(
                f"barycorr infile={i} outfile=bc{evttype}_{obsid}_{ccn}.evt"
                + f" orbitfiles={orbfile} refframe=ICRS"
                + f" ra={self.ras[index]} dec={self.decs[index]}",
                shell=True,
            )
            comp_org = compress_gz(i)
            print(comp_org)
            ccn = ccn + 1
        os.chdir(self.base_dir)

    def pull_reduce(self, bc):
        """
        Pulls all selected OBSIDs
        """
        for i in self.obs:
            index = self.obs.index(i)
            self.pull(index)
            if bc is True:
                self.barycenter(index)


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
    args = cli()
    xte = Autoxte(target=args.source)
    xte.selection()
    xte.pull_reduce(args.bc)
