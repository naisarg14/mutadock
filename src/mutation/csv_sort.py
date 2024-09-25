################################################################################
#                           PROJECT INFORMATION                                #
#                              Name: MUTADOCK                                  #
#                           Author: Naisarg Patel                              #
#                                                                              #
#       Copyright (C) 2024 Naisarg Patel (https://github.com/naisarg14)        #
#                                                                              #
#          Project: https://github.com/naisarg14/mutadock                      #
#                                                                              #
#   This program is free software; you can redistribute it and/or modify it    #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
################################################################################


from mutation.helpers import backup
import sys, os
import argparse

try:
    import pandas as pd
except ImportError:
    msg = "Error with importing pandas module for mutation using mutadock.\n"
    msg += "Easiest way to fix this is to install pandas using the following command:\n\n"
    msg += "python -m pip install pandas\n"
    msg += "If you already have pandas installed, please check the installation.\n"
    msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
    print(msg)
    sys.exit(2)


def main():
    parser = argparse.ArgumentParser(description="This program takes as input a CSV file and sorts it according to the coloumn given.", epilog="Written by Naisarg Patel (https://github.com/naisarg14)")
    parser.add_argument("-i",'--input', help="CSV File to sort", metavar="CSV", required=True)
    parser.add_argument("-o", "--output", help="Sorted Output CSV", metavar="CSV")

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-a', action='store_true', help='Sort in ascending order')
    group.add_argument('-d', action='store_true', help='Sort in descending order')

    parser.add_argument("-n", "--name", metavar="NAME", help="Name of the coloumn to sort")
    parser.add_argument("-N", "--num", metavar="NUM", help="Number of the coloumn to sort")

    args = parser.parse_args()
    if args.d:
        order = False
    else:
        order = True

    in_file = args.input


    if not os.path.isabs(in_file):
        print("Assuming current directory as root since path not specified.")
        dir = os.getcwd()
        file = in_file
    else:
        dir, file = os.path.split(os.path.abspath(in_file))

    full_csv_path = os.path.join(dir, file)

    if not os.path.isfile(full_csv_path):
        sys.exit("No such file found in current directory, enter full path for other directories.")

    if not file.endswith(".csv"):
        sys.exit("Given file is not a CSV file, input should be a CSV file.")

    sort_csv(in_file=full_csv_path, out_file=args.output, col_num=args.num, col_name=args.name, order=order)


def sort_csv(in_file, out_file=None, col_num=None, col_name=None, order=True):
    in_file = in_file.removesuffix('.csv')
    if not out_file:
        out_file = f"{in_file.removesuffix('.csv')}_sorted.csv"
    backup(out_file)
    try:
        df = pd.read_csv(f"{in_file}.csv")
    except FileNotFoundError:
        print("No such file found in current directory, enter full path for other directories.")
        return 420

    if col_name:
        header = list(df.columns.values)
        for i in range(len(header)):
            if header[i].lower() == col_name.lower():
                col_num = i
                break

    if not col_num:
        count = 0
        for i in (df.columns.values):
            print(f"{count}: {i}")
            count += 1
        col_num = int(input("Enter the column number to sort by: "))
    df = df.sort_values(by=df.columns[col_num], ascending=True)
    df = df.iloc[: , 1:]
    df.insert(0, 'sr', range(1, 1 + df.shape[0]))
    df.to_csv(out_file, index=False)

    return out_file



if __name__ == "__main__":
    main()
