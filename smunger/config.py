# Columns
FORMATTED_COLS = ["#CHR", "BP", "rsID", "EA", "NEA", "EAF", "MAF", "BETA", "SE", "P"]

# Data types
COL_TYPES = {
    "#CHR": int,
    "BP": int,
    "rsID": str,
    "EA": str,
    "NEA": str,
    "EAF": float,
    "MAF": float,
    "BETA": float,
    "SE": float,
    "P": float,
}

# allow for missing values
COL_ALLOWNA = {
    "#CHR": False,
    "BP": False,
    "rsID": True,
    "EA": False,
    "NEA": False,
    "EAF": True,
    "MAF": True,
    "BETA": True,
    "SE": True,
    "P": True,
}

# data ranges
COL_RANGE = {
    "#CHR": (1, 22),
    "BP": (0, 999999999),
    "rsID": None,
    "EA": None,
    "NEA": None,
    "EAF": (0, 1),
    "MAF": (0, 1),
    "BETA": None,
    "SE": (0, None),
    "P": (0, 1),
}
