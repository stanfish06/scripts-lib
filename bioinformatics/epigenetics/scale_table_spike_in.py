#!/usr/bin/env python
import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--fname", required=True)
parser.add_argument("--spike_scale", type=float, default=100000)
args = parser.parse_args()
df = pd.read_csv(args.fname)
df["fraction"] = args.spike_scale * df["n_mapped"] / df["total"]
df["scale_factor"] = 1 / df["fraction"]
df.to_csv("scale_table.csv", index=False)
