import os
import pandas as pd

path = '/disk/raptor/lclhome/mtari008/colab/SORI_imput_files'
files = os.listdir(path)

files_xls = [os.path.join(path, f) for f in files if f[-5:] == '.xlsx']

df = pd.DataFrame()

for f in files_xls:
    data = pd.read_excel(f, 'Sheet1')
    df = df.append(data)

df.reset_index(drop=True, inplace=True)
# pd.set_option("display.max_rows", None, "display.max_columns", None)
print(df)
df.to_excel("input/SORI-combined.xlsx", index=False)