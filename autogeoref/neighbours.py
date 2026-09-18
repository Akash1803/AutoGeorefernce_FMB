"""Which survey numbers each FMB sheet prints outside its own outline, and on which side.

Source of truth is the two-reader transcription at
FMB_Vector/<village>/neighbour_transcription/nb_override_<village>.json. The glyph reader of
2026-09-17 misread 170 as "17" and took the adjoining village label "V.No. 74 Thirukachur" for a
survey, so it is not used here.
"""
import json
import re

from . import paths

SIDE_VEC = {"N": (0.0, 1.0), "NE": (0.7071, 0.7071), "E": (1.0, 0.0), "SE": (0.7071, -0.7071),
            "S": (0.0, -1.0), "SW": (-0.7071, -0.7071), "W": (-1.0, 0.0), "NW": (-0.7071, 0.7071)}
_CACHE = {}


def load(village):
    if village not in _CACHE:
        p = paths.vector_dir(village) / "neighbour_transcription" / ("nb_override_%s.json" % village)
        _CACHE[village] = json.loads(p.read_text(encoding="utf-8")) if p.exists() else None
    return _CACHE[village]


def normalise(number):
    """'103/5A' -> '103'; 'V.No. 74 ...' -> None (an adjoining village, not a survey)."""
    s = str(number).upper().strip()
    if "V.NO" in s or "VILLAGE" in s:
        return None
    s = re.sub(r"^S\.?\s*NO\.?\s*", "", s).replace(" ", "").split("/")[0]
    return s if re.fullmatch(r"\d+[A-Z]?", s) else None


def printed(village, survey):
    """Survey ids this sheet prints outside its outline, trusted entries only."""
    table = load(village)
    if not table or survey not in table:
        return set()
    out = set()
    for entry in table[survey]:
        n = normalise(entry.get("number", ""))
        if not n or n == str(survey).upper():
            continue
        if entry.get("readers", 1) >= 2 or entry.get("confidence") == "high":
            out.add(n)
    return out


def are_neighbours(village, a, b):
    """True / False / None (unknown, because at least one sheet was not transcribed)."""
    table = load(village)
    if not table or a not in table or b not in table:
        return None
    return str(b).upper() in printed(village, a) or str(a).upper() in printed(village, b)


def side(village, survey, other):
    table = load(village)
    if not table or survey not in table:
        return None
    for entry in table[survey]:
        if normalise(entry.get("number", "")) == str(other).upper():
            s = str(entry.get("side", "")).upper()
            return s if s in SIDE_VEC else None
    return None


def side_vector(letters):
    return SIDE_VEC.get(str(letters).upper(), (0.0, 0.0))
