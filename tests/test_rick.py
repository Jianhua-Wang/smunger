import pandas as pd
from rich.table import Table
from smunger.console import console
from smunger.smunger import map_colnames

df = pd.read_csv('./exampledata/test.txt', sep='\t')
df = df.rename(columns={'CHR': 'chr1', 'BP': 'pos'})
# a = df.copy()
# a.columns = [str(i) for i in range(len(a.columns))]
# b = df.copy()
# b.columns = [str(i)+'a' for i in range(len(b.columns))]
# df = pd.concat([df, a, b], axis=1)

colname_map = {
    'CHR': 'CHR',
    'BETA': 'BETA',
}
if len(df) < 5:
    nrows = len(df)
else:
    nrows = 5

map_colnames(df)
