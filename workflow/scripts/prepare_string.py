import argparse
import pandas as pd
import sys



parser = argparse.ArgumentParser()

parser.add_argument('--input', type=str, help='STRING data', required=True)
parser.add_argument('--output', type=str, help='Processed STRING network', required=True)


args = parser.parse_args()

print("Python version:")
print(sys.version)
print("Arguments:")
print(args)
print("\n")


input_file = args.input
output_file = args.output

df = pd.read_csv(input_file)

print("Unprocessed STRING data:")
print(df)


df = df[["protein1","protein2"]]

# remove duplicates
df  = df.drop_duplicates()

print("\nSTRING network without duplicates:")
print(df)

# remove self-regulating edges

df = df.loc[df["protein1"] != df["protein2"]]

print("\nProcessed STRING network without self-regulations:")
print(df)

df.to_csv(output_file, index=False)