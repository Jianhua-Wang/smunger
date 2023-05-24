"""Define the constants used in the package."""

import numpy as np


class ColName:
    """Define column names."""

    CHR = "CHR"
    BP = "BP"
    RSID = "rsID"
    SNPID = 'SNPID'  # unique snpid, chr-bp-sorted(EA,NEA)
    EA = "EA"
    NEA = "NEA"
    P = "P"
    NEGLOGP = "-log10P"
    BETA = "BETA"
    OR = "OR"
    ORSE = "ORSE"
    SE = "SE"
    EAF = "EAF"
    MAF = "MAF"
    N = "N"
    Z = "Z"
    INFO = "INFO"
    OUTCOLS = ['CHR', 'BP', 'rsID', 'EA', 'NEA', 'EAF', 'MAF', 'BETA', 'SE', 'P']


chrom_len = {
    1: 249250621,
    2: 243199373,
    3: 198022430,
    4: 191154276,
    5: 180915260,
    6: 171115067,
    7: 159138663,
    8: 146364022,
    9: 141213431,
    10: 135534747,
    11: 135006516,
    12: 133851895,
    13: 115169878,
    14: 107349540,
    15: 102531392,
    16: 90354753,
    17: 81195210,
    18: 78077248,
    19: 59128983,
    20: 63025520,
    21: 48129895,
    22: 51304566,
}


class ColType:
    """Define column types."""

    CHR = int
    BP = int
    RSID = str
    EA = str
    NEA = str
    P = float
    NEGLOGP = float
    BETA = float
    OR = float
    ORSE = float
    SE = float
    EAF = float
    MAF = float
    N = int
    Z = float
    INFO = float


class ColAllowNA:
    """Define whether a column allows missing values."""

    CHR = False
    BP = False
    RSID = True
    EA = False
    NEA = False
    P = True
    BETA = True
    OR = True
    ORSE = True
    SE = True
    EAF = True
    MAF = True
    N = True
    Z = True
    INFO = True


class ColRange:
    """Define the range of values for each column."""

    INFO = (0, 1)
    CHR_MIN = 1
    CHR_MAX = 23
    BP_MIN = 0
    BP_MAX = 300000000
    P_MIN = 0
    P_MAX = 1
    SE_MIN = 0
    SE_MAX = np.inf
    EAF_MIN = 0
    EAF_MAX = 1
    MAF_MIN = 0
    MAF_MAX = 1
    ORSE_MIN = 0
    ORSE_MAX = np.inf
    NEGLOGP_MIN = 0


COMMON_COLNAMES = {
    'CHR': 'CHR',
    'chr': 'CHR',
    'chromosome': 'CHR',
    '#CHROM': 'CHR',
    '#CHR': 'CHR',
    'Chromosome': 'CHR',
    'CHROM': 'CHR',
    'Chr': 'CHR',
    'BP': 'BP',
    'position': 'BP',
    'pos': 'BP',
    'bp': 'BP',
    'base_pair_location': 'BP',
    'POS': 'BP',
    'Position': 'BP',
    'Position_b37': 'BP',
    'Bp': 'BP',
    'MarkerName': 'rsID',
    'varID': 'rsID',
    'rsid': 'rsID',
    'snp': 'rsID',
    'SNP': 'rsID',
    'hm_rsid': 'rsID',
    'variant_id': 'rsID',
    'uniqID': 'rsID',
    'ID': 'rsID',
    'id': 'rsID',
    'rsID': 'rsID',
    'Variant_ID': 'rsID',
    'Test': 'rsID',
    '#SNP': 'rsID',
    'MAF': 'MAF',
    'all_maf': 'MAF',
    'beta': 'BETA',
    'BETA_INSOMNIA': 'BETA',
    'BETA': 'BETA',
    'all_inv_var_meta_beta': 'BETA',
    'Effect': 'BETA',
    'lnOR': 'BETA',
    'Beta': 'BETA',
    'LogOR': 'BETA',
    'frequentist_add_beta_1': 'BETA',
    'se': 'SE',
    'or_se': 'SE',
    'standard_error': 'SE',
    'SE_INSOMNIA': 'SE',
    'SE': 'SE',
    'sebeta': 'SE',
    'all_inv_var_meta_sebeta': 'SE',
    'StdErr': 'SE',
    'Standard_error': 'SE',
    'SE_GC': 'SE',
    'StdErrLogOR': 'SE',
    'frequentist_add_se_1': 'SE',
    'P-value': 'P',
    'Score.pval': 'P',
    'p': 'P',
    'p.metal': 'P',
    'P': 'P',
    'p_value': 'P',
    'P_INSOMNIA': 'P',
    'P_BOLT_LMM_INF': 'P',
    'p.value': 'P',
    'pval': 'P',
    'all_inv_var_meta_p': 'P',
    'PS': 'P',
    'P_Value': 'P',
    'Pvalue': 'P',
    'P_GC': 'P',
    'P.value': 'P',
    'Zscore': 'Zscore',
    'Z': 'Zscore',
    'z_value': 'Zscore',
    'z.meta': 'Zscore',
    'Z-score': 'Zscore',
    'Freq1': 'EAF',
    'altFreq1': 'EAF',
    'af': 'EAF',
    'freqa1': 'EAF',
    'effect_allele_frequency': 'EAF',
    'A1FREQ': 'EAF',
    'all_meta_AF': 'EAF',
    'FREQ_Effect_Allele': 'EAF',
    'EAF': 'EAF',
    'Frq': 'EAF',
    'Coded_freq': 'EAF',
    'Freq': 'EAF',
    'eaf': 'EAF',
    'or': 'OR',
    'odds_ratio': 'OR',
    'OR': 'OR',
    'ORX': 'OR',
    'ORS': 'OR',
    'Odds_ratio': 'OR',
    'other_allele': 'NEA',
    'effect_allele': 'EA',
    'Allele1': 'EA',
    'Allele2': 'NEA',
    'Freq.Allele1.HapMapCEU': 'EAF',
    'allele1': 'EA',
    'allele2': 'NEA',
    'freqA1': 'EAF',
    'A1': 'EA',
    'A2': 'NEA',
    'maf': 'MAF',
    'beta_SNP_add': 'BETA',
    'sebeta_SNP_add': 'SE',
    'Minor_allele': 'EA',
    'all_OR': 'OR',
    'frequentist_add_pvalue': 'P',
    'z': 'Zscore',
    'SNP_ID': 'rsID',
    'Allele 1': 'EA',
    'Allele 2': 'NEA',
    'Effect allele (EA)': 'EA',
    'Effect allele frequency (EAF)': 'EAF',
    'P value': 'P',
    'AF': 'EAF',
    'ALLELE1': 'EA',
    'ALLELE0': 'NEA',
    'P_LINREG': 'P',
    'snpid': 'rsID',
    'bpos': 'BP',
    'a1': 'EA',
    'a2': 'NEA',
    'mtag_beta': 'BETA',
    'mtag_se': 'SE',
    'mtag_z': 'Zscore',
    'mtag_pval': 'P',
    'position_GRCh37': 'BP',
    'non_effect_allele': 'NEA',
    'REF': 'NEA',
    'ALT': 'EA',
    'noneffect_allele': 'NEA',
    'FreqSE': 'SE',
    'ALT_ALLELE_FREQ': 'EAF',
    'PVALUE': 'Zscore',
    'EAF_A1': 'EAF',
    'Pval': 'P',
    'RSID': 'rsID',
    'a0': 'NEA',
    'freq1': 'EAF',
    'freq_se': 'SE',
    'EA': 'EA',
    'NEA': 'NEA',
    'LOGOR': 'BETA',
    'LOGOR_SE': 'SE',
    'Effect_allele': 'EA',
    'Non_Effect_allele': 'NEA',
    'PVAL': 'P',
    'Effect_Allele': 'EA',
    'NonEffect_Allele': 'NEA',
    'Effect_Allele_Frequency': 'EAF',
    'SE_of_Beta': 'SE',
    'A1Frq': 'EAF',
    'FREQ': 'EAF',
    'a_0': 'NEA',
    'a_1': 'EA',
    'pvalue': 'P',
    'Pos_b37': 'BP',
    'Allele_1': 'EA',
    'Allele_2': 'NEA',
    'Other_allele': 'NEA',
    'Effect_allele_Freq': 'EAF',
    'markername': 'rsID',
    'effect_allele_freq': 'EAF',
    'EFFECT': 'BETA',
    'STDERR': 'SE',
    'OA': 'NEA',
    'BPos': 'BP',
    'Non_effect_Allele': 'NEA',
    'Standard_Error_of_Beta': 'SE',
    'p-value': 'P',
    'snp_ids': 'rsID',
    'Ref_allele': 'NEA',
    'Alt_allele': 'EA',
    'JASS_PVAL': 'P',
    'variant_ID': 'rsID',
    'alleleA': 'EA',
    'alleleB': 'NEA',
    'A1Freq': 'EAF',
    'MAF_calculated_from_dosage_data': 'MAF',
    'REF_allele': 'NEA',
    'ALT_allele': 'EA',
    'none_effect_alllele': 'NEA',
    'Pos': 'BP',
    'ZSCORE': 'Zscore',
    'zscore': 'Zscore',
    'allele_1': 'EA',
    'allele_2': 'NEA',
    'allele0': 'NEA',
    'snp.1': 'rsID',
    'a1freq': 'EAF',
    'snpname': 'rsID',
    'ref': 'NEA',
    'alt': 'EA',
    'A0': 'NEA',
    'reference_allele': 'NEA',
    'base_pair_location_grch37': 'BP',
    'Rsid': 'rsID',
    'Non_effect_allele': 'NEA',
    'EffectAllele': 'EA',
    'NonEffectAllele': 'NEA',
    'Freq.A1.ESP.EUR': 'EAF',
    'base_pair_position': 'BP',
    'effect_allle_frequency': 'EAF',
    'base_pair_locations': 'BP',
    'ZScore': 'Zscore',
    'stderr': 'SE',
    'pval(-log10)': '-log10P',
    '-log10P': '-log10P',
    '-log10p': '-log10P',
    '-log10Pvalue': '-log10P',
}
