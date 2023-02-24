import tabix
import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

file = '/f/jianhua/Resources/rs2pos/pos2rsid.txt.gz'
tb = tabix.open(file)
outfile = '~/catalog.munged.rsid.txt'

chunksize = 2000000
ith = 0
for df in pd.read_csv('~/catalog.munged.txt.gz', sep='\t', chunksize=100000):
    for chrom, chr_df in df.groupby('CHR'):
        chr_df = chr_df.sort_values('BP')
        for i in range(chr_df['BP'].min(), chr_df['BP'].max(), chunksize):
            chunk_df = chr_df[(chr_df['BP'] >= i) & (chr_df['BP'] < i + chunksize)].copy()
            if len(chunk_df) == 0:
                continue
            rsid_map = pd.DataFrame(
                columns=['chr', 'bp', 'rsid', 'ref', 'alt'], data=tb.query('1', i - 1, i + chunksize)
            )
            rsid_map = rsid_map.drop_duplicates(subset=['bp'])
            rsid_map['bp'] = rsid_map['bp'].astype(int)
            rsid_map = pd.Series(data=rsid_map['rsid'].values, index=rsid_map['bp'].values)
            chunk_df['rsID'] = chunk_df['BP'].map(rsid_map)
            logging.info(f'Processing {chrom}:{i}-{i + chunksize}, chunk No.{ith}, {len(chunk_df)} rows.')
            ith += 1
            if ith == 1:
                chunk_df.to_csv(outfile, sep='\t', index=False, mode='w')
            else:
                chunk_df.to_csv(outfile, sep='\t', index=False, mode='a', header=False)
