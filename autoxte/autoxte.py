#! /usr/bin/env python3

# AutoXTE
# Copyright 2023 Nicholas Kuechel
# License Apache 2.0


import subprocess as sp
import os
import argparse
from astroquery.heasarc import Heasarc
from astropy.table import Table
from astropy.time import Time
import wget
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


class Autoxte:
    def __init__(self):
        self.state = True
        self.obs = []
        self.ras = []
        self.decs = []
        self.cNums = []
        self.prnbs = []
        self.base_dir = os.getcwd()
        print("##############  Auto XTE  ##############")

    def query(self, target=None):
        """
        Querys xtemaster @ heasarc for a specific target
        """
        obj = str(input("Target: "))
        heasarc = Heasarc()
        pca = heasarc.query_object(obj, mission="xtemaster", resultmax=1000000)
        pca = Table(pca).to_pandas()
        cnt = 0
        for i in pca["OBSID"]:
            i = i.decode()
            pca.loc[cnt, "OBSID"] = str(i)
            cnt = cnt + 1
        cnt = 0
        for i in pca["TIME"]:  # converts times from mjd to datetime format
            t0 = Time(i, format="mjd").to_datetime()
            pca.loc[cnt, "TIME"] = t0
            cnt = cnt + 1
        pca = pca.loc[pca["EXPOSURE"] != 0]
        return pca

    def selection(self, table):
        """
        Interface for selecting the OBSIDs to be retrieved
        """
        while self.state is True:
            enter = str(input("autoXTE >> "))
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
                    rows = table.loc[table["PRNB"] == int(enter)]
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
                    for i in table["OBSID"]:
                        self.obs.append(i)
                    for i in table["RA"]:
                        self.ras.append(i)
                    for i in table["DEC"]:
                        self.decs.append(i)
                    for i in table["PRNB"]:
                        self.prnbs.append(i)
                        prnb = str(i)
                        if int(prnb[0]) == 9:
                            cycles = [9, 10, 11, 12, 13, 14, 15, 16]
                            self.cNums.append(cycles[int(prnb[1])])
                        else:
                            self.cNums.append(int(prnb[0]))
                else:
                    self.obs.append(enter)
                    row = table.loc[table["OBSID"] == enter]
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

    def retrieve(self, index):
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
        wget.download(fullurl)

    def barycenter(self, index):
        """
        Calls the HEASoft's barycorr tool to perform
        a barycenter correction to .evt data
        """
        obsid = self.obs[index]
        orb_file = glob.glob(f"P{self.prnbs[index]}/" +
                             f"{self.obs[index]}/orbit/FP*")
        orb_file = f"{self.base_dir}/{orb_file[0]}"
        os.chdir("P{self.prnbs[index]}/{self.obs[index]}/pca")
        gzevts = glob.glob("*.evt.gz")
        for i in gzevts:
            extract_gz(i)
        evt_files = glob.glob("*evt")
        ccn = 1
        for i in evt_files:
            sp.call(
                f"barycorr infile={i} outfile=bcSE_{obsid}_{ccn}.evt"
                + f"orbitfiles={orb_file} refframe=ICRS"
                + f"ra={self.ras[index]} dec={self.decs[index]}",
                shell=True,
            )
            ccn = ccn + 1
        os.chdir(self.base_dir)


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
