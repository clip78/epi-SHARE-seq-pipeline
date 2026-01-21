#!/usr/bin/env python3

"""
This script takes in the Picard CollectInsertSizeMetrics histogram txt file output,
and generates the histogram as a png.
"""

import argparse
import pandas as pd
from plotnine import *

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot insert size histogram")
    parser.add_argument("histogram_file", help="Histogram txt file name")
    parser.add_argument("pkr", help="PKR ID")
    parser.add_argument("out_file", help="Name of output png file")
    
    return parser.parse_args()

def get_hist_vals(histogram_file):
    """Get dataframe of histogram values"""
    df = pd.read_csv(histogram_file, sep="\t", comment="#")
    df = df.rename(columns={df.columns[1]: "count"})
    df = df[["insert_size", "count"]]
    return(df)

def label_func(breaks):
    return ["{:.0e}".format(x) for x in breaks]

def plot_hist(df, pkr, out_file):
    plot = (ggplot(df, aes(x="insert_size", y="count")) +
            geom_line(color="red") +
            geom_area(fill="red") +
            labs(title = f"Insert Size Histogram ({pkr})",
                 x = "Insert size",
                 y = "Count") + 
            scale_y_continuous(labels = label_func) +
            theme_classic())
    
    plot.save(filename = out_file, dpi=1000)

def main():
    print("Starting histogram plotting script")
    args = parse_arguments() 
    histogram_file = getattr(args, "histogram_file")
    pkr = getattr(args, "pkr")
    out_file = getattr(args, "out_file")
    
    df = get_hist_vals(histogram_file)
    
    plot_hist(df, pkr, out_file)
    print("Finished plotting")

if __name__ == "__main__":
    main()    
