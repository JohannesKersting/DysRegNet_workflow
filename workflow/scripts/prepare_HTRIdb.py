import argparse
import pandas as pd
import sys


parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str, help='Unprocessed HTRIdb data', required=True)
parser.add_argument('--output', type=str, help='Processed HTRIdb network', required=True)


args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")


input_file = args.input
output_file = args.output

df = pd.read_csv(input_file, sep=";")

print("Unprocessed HTRIdb data:")
print(df)


df = df[["SYMBOL_TF","SYMBOL_TG"]]

# remove duplicates
df  = df.drop_duplicates()

print("\nHTRIdb network without duplicates:")
print(df)

# remove self-regulating edges

df = df.loc[df["SYMBOL_TF"] != df["SYMBOL_TG"]]

print("\nProcessed HTRIdb network without self-regulations:")
print(df)

df.to_csv(output_file, index=False)