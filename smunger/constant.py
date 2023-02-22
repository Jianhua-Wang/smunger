"""Define the constants used in the package."""


class ColName:
    """Define column names."""

    CHR = "CHR"
    BP = "BP"
    RSID = "rsID"
    EA = "EA"
    NEA = "NEA"
    P = "P"
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


class ColType:
    """Define column types."""

    CHR = int
    BP = int
    RSID = str
    EA = str
    NEA = str
    P = float
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

    CHR = (0, 24)
    BP = (0, 300000000)
    RSID = None
    EA = None
    NEA = None
    P = (0, 1)
    BETA = None
    OR = None
    ORSE = None
    SE = (0, None)
    EAF = (0, 1)
    MAF = (0, 1)
    N = (0, None)
    Z = None
    INFO = (0, 1)
