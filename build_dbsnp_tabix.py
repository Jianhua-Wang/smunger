"""Build a tabix-indexed database of dbSNP data."""
import json
import bz2
import numpy as np
import pandas as pd
import gzip
from subprocess import call, check_output


version = 'b156'

# Download the dbSNP data
# call(f'wget https://ftp.ncbi.nlm.nih.gov/snp/archive/{version}/VCF/GCF_000001405.25.gz', shell=False)
# call(f'wget https://ftp.ncbi.nlm.nih.gov/snp/archive/{version}/VCF/GCF_000001405.25.gz.tbi', shell=False)
# call(f'wget ftp://ftp.ncbi.nlm.nih.gov/snp/archive/{version}/JSON/refsnp-merged.json.bz2', shell=True)


# split multi-allelic variants into separate rows
call('bcftools norm -m - GCF_000001405.25.gz -O z -o GCF_000001405.25.nomultiallelic.gz', shell=True)


n_header = 0
with gzip.open('./GCF_000001405.25.nomultiallelic.gz', 'rb') as f:
    for line in f:
        if line.decode('utf-8').startswith('##'):
            n_header += 1
        else:
            break

chr_name_map = pd.Series(
    data=[str(i) for i in range(1, 23)] + ['X', 'Y'],
    index=check_output('tabix -l ./GCF_000001405.25.gz', shell=True).decode().split()[:24],
)

for df in pd.read_csv(
    './GCF_000001405.25.nomultiallelic.gz',
    sep='\t',
    skiprows=n_header,
    compression='gzip',
    usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT'],
    chunksize=100000,
):
    # replace chr
    df['#CHROM'] = df['#CHROM'].map(chr_name_map)

    # remove chrM
    df = df.dropna()
    # stop when df are all chrM
    if len(df) == 0:
        break

    # use first int of rsid as fake chr
    df['rsid'] = df['ID'].map(lambda rsid: rsid[2:])
    df['rsid_1st'] = df['ID'].map(lambda rsid: rsid[2])

    # write pos2snp file
    df[['#CHROM', 'POS', 'ID', 'REF', 'ALT']].to_csv(
        f'./pos2snp_{version}.txt', sep='\t', index=False, header=False, mode='a'
    )
    # write snp2pos file
    df[['rsid_1st', 'rsid', '#CHROM', 'POS', 'REF', 'ALT']].to_csv(
        f'./snp2pos_{version}_unsorted.txt', sep='\t', index=False, header=False, mode='a'
    )


# sort snp2pos file
call(f'sort -k1,1 -k2,2n ./snp2pos_{version}_unsorted.txt > ./snp2pos_{version}.txt', shell=True)

# bgzip and tabix
call(f'bgzip ./pos2snp_{version}.txt', shell=True)
call(f'tabix -s 1 -b 2 -e 2 ./pos2snp_{version}.txt.gz', shell=True)
call(f'bgzip ./snp2pos_{version}.txt', shell=True)
call(f'tabix -s 1 -b 2 -e 2 ./snp2pos_{version}.txt.gz', shell=True)

# parse merged rsids
with open(f'./merged_{version}.txt', 'w') as f_out:
    with bz2.BZ2File('./refsnp-merged.json.bz2', 'rb') as f_in:
        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            refsnp_id = rs_obj['refsnp_id']
            for merge_into in rs_obj['merged_snapshot_data']['merged_into']:
                f_out.write(f"{refsnp_id[0]}\t{refsnp_id}\t{merge_into}\n")
                for dbsnp1_merges in rs_obj['dbsnp1_merges']:
                    dbsnp1_merges = dbsnp1_merges['merged_rsid']
                    f_out.write(f"{dbsnp1_merges[0]}\t{dbsnp1_merges}\t{merge_into}\n")


# bgzip and tabix
call(f'bgzip ./merged_{version}.txt', shell=True)
call(f'tabix -s 1 -b 2 -e 2 ./merged_{version}.txt.gz', shell=True)

# remove intermediate files
call('rm ./GCF_000001405.25.gz', shell=True)
call('rm ./GCF_000001405.25.nomultiallelic.gz', shell=True)
call('rm ./refsnp-merged.json.bz2', shell=True)
call(f'rm ./snp2pos_{version}_unsorted.txt', shell=True)
