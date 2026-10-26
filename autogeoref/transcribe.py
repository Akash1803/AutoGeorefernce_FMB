"""Two-reader transcription of the survey numbers printed around each FMB sheet (PR 1).

The engine trusts a printed neighbour only when two independent readers agree on it
(``neighbours.printed``). This module does everything around those readers: it renders the
drawing area of every sheet to a PNG, loads two finished readings, merges them, writes the
engine's table ``nb_override_<village>.json``, lists every disagreement in a review CSV for the
analyst, writes an agreement report with a random sample of agreed sheets to spot-check, and
computes the seed set (greedy dominating set) over the resulting neighbour graph.

Reader backends are named in the config (``transcription.readers``, exactly two). ``json:<name>``
loads a finished reading from ``<vector_dir>/neighbour_transcription/readings/<name>.json``,
written by a vision-model session that looked at the rendered crops. No OCR runs here and
nothing is trained; PaddleOCR is PR 4 and will be another backend name.

Both readers in use on 2026-09-21 are the same vision model in separate contexts, so their
agreement overstates accuracy against the PDF; the analyst's sample check is the accuracy figure.
"""
import csv
import dataclasses
import datetime
import json
import logging
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, List, Optional, Sequence, Set, Tuple

from . import config as configmod, neighbours, paths

log = logging.getLogger(__name__)

SIDES = tuple(neighbours.SIDE_VEC)
REVIEW_COLUMNS = ["village", "survey", "field", "number", "reader_a", "printed_a", "side_a",
                  "reader_b", "printed_b", "side_b", "resolution", "resolved_by", "note"]
SOURCE_READER, SOURCE_HAND = "reader", "hand"
RESOLVE_YES = ("yes", "y", "keep", "ok", "true", "1")
RESOLVE_NO = ("no", "n", "drop", "remove", "false", "0")


@dataclass(frozen=True)
class Entry:
    """One printed label as a reader saw it."""
    number: str                      # as printed, e.g. "569/C"
    side: str = ""                   # compass side of the label around the outline, or ""
    confidence: str = "high"         # high | low


Reading = Dict[str, List[Entry]]     # survey -> labels read on that sheet


@dataclass
class Disagreement:
    village: str
    survey: str
    field: str                       # "number" (seen by one reader) or "side"
    number: str                      # normalised survey id
    reader_a: str
    printed_a: str
    side_a: str
    reader_b: str
    printed_b: str
    side_b: str
    resolution: str = ""
    resolved_by: str = ""            # "hand" once the analyst filled the resolution
    note: str = ""

    def key(self) -> Tuple[str, str, str]:
        return (self.survey, self.field, self.number)


@dataclass
class SheetStats:
    survey: str
    n_a: int
    n_b: int
    n_both: int
    n_union: int
    sides_checked: int
    sides_agreed: int

    @property
    def jaccard(self) -> Optional[float]:
        return round(self.n_both / self.n_union, 3) if self.n_union else None


# ----------------------------------------------------------------------------- rendering


def render_dir(cfg: Optional[configmod.Config] = None) -> Path:
    """Where sheet crops go: gitignored cache under the project, deletable at any time."""
    cfg = cfg or configmod.load()
    custom = cfg.transcription.render_dir
    return Path(custom) if custom else paths.PROJECT / "_cache" / "sheet_renders"


def render_sheet(pdf_path: Path, out_path: Path, scale: float = 1.6, pad_px: int = 120,
                 ink_below: int = 160) -> Path:
    """Render page 1 in grey and crop to the drawing: header and footer bands are skipped, the
    bounding box of dark pixels is padded. Returns ``out_path`` (unchanged if already newer)."""
    import numpy as np
    import pymupdf
    from PIL import Image
    pdf_path, out_path = Path(pdf_path), Path(out_path)
    if out_path.exists() and out_path.stat().st_mtime >= pdf_path.stat().st_mtime:
        return out_path
    doc = pymupdf.open(str(pdf_path))
    page = doc[0]
    pix = page.get_pixmap(matrix=pymupdf.Matrix(scale, scale), colorspace=pymupdf.csGRAY)
    a = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width)
    h, w = a.shape
    y_lo, y_hi, x_lo, x_hi = int(h * 0.05), int(h * 0.95), int(w * 0.02), int(w * 0.98)
    body = a[y_lo:y_hi, x_lo:x_hi]
    ys, xs = np.where(body < ink_below)
    if len(ys) == 0:
        crop = a
    else:
        y0, y1 = max(0, ys.min() + y_lo - pad_px), min(h, ys.max() + y_lo + pad_px)
        x0, x1 = max(0, xs.min() + x_lo - pad_px), min(w, xs.max() + x_lo + pad_px)
        crop = a[y0:y1, x0:x1]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(crop).save(str(out_path))
    return out_path


def render_village(village: str, cfg: Optional[configmod.Config] = None) -> List[Path]:
    """Crops for every sheet of the village, under ``render_dir(cfg)/<village>/<survey>.png``."""
    cfg = cfg or configmod.load()
    out_dir = render_dir(cfg) / village
    done = []
    for survey in paths.surveys_with_sheets(village):
        pdf = paths.sketch_dir(village) / ("%s.pdf" % survey)
        if not pdf.exists():
            log.warning("%s %s: no PDF at %s", village, survey, pdf)
            continue
        done.append(render_sheet(pdf, out_dir / ("%s.png" % survey), scale=cfg.transcription.render_scale))
    log.info("%s: %d sheet crops under %s", village, len(done), out_dir)
    return done


# ----------------------------------------------------------------------------- readers


def readings_dir(village: str) -> Path:
    return paths.vector_dir(village) / "neighbour_transcription" / "readings"


class JsonReader:
    """A finished reading on disk: ``{survey: [{"number", "side", "confidence"}, ...]}``."""

    def __init__(self, name: str):
        self.name = name

    def path(self, village: str) -> Path:
        return readings_dir(village) / ("%s.json" % self.name)

    def read(self, village: str) -> Reading:
        p = self.path(village)
        if not p.exists():
            raise FileNotFoundError("reader %s has no reading for %s at %s" % (self.name, village, p))
        return parse_reading(json.loads(p.read_text(encoding="utf-8")), where=str(p))


def parse_reading(data: Dict[str, List[dict]], where: str = "") -> Reading:
    if not isinstance(data, dict):
        raise ValueError("reading %s must be an object keyed by survey" % where)
    out: Reading = {}
    for survey, entries in data.items():
        if not isinstance(entries, list):
            raise ValueError("reading %s: %s must hold a list" % (where, survey))
        items = []
        for e in entries:
            number = str(e.get("number", "")).strip()
            if not number:
                continue
            side = str(e.get("side", "") or "").upper().strip()
            if side and side not in SIDES:
                raise ValueError("reading %s: %s prints %r with side %r; expected one of %s" % (where, survey, number, side, SIDES))
            conf = str(e.get("confidence", "high") or "high").lower()
            items.append(Entry(number, side, "low" if conf == "low" else "high"))
        out[str(survey)] = items
    return out


def reader_for(spec: str) -> JsonReader:
    kind, _, name = spec.partition(":")
    if kind != "json" or not name:
        raise ValueError("reader %r: only 'json:<name>' backends exist in PR 1" % spec)
    return JsonReader(name)


def load_readings(village: str, cfg: Optional[configmod.Config] = None) -> Dict[str, Reading]:
    cfg = cfg or configmod.load()
    return {spec: reader_for(spec).read(village) for spec in cfg.transcription.readers}


# ----------------------------------------------------------------------------- merge


def _by_number(entries: Iterable[Entry]) -> Dict[str, List[Entry]]:
    out: Dict[str, List[Entry]] = {}
    for e in entries:
        n = neighbours.normalise(e.number)
        if n is None:
            continue                                  # a village label or noise, not a survey
        out.setdefault(n, []).append(e)
    return out


def merge(reading_a: Reading, reading_b: Reading, surveys: Optional[Iterable[str]] = None,
          village: str = "", name_a: str = "A", name_b: str = "B"
          ) -> Tuple[Dict[str, List[dict]], List[Disagreement], List[SheetStats]]:
    """Two readings -> (engine table, disagreements, per-sheet agreement).

    A number both readers print is trusted (``readers: 2``); its side is the side both gave, or
    empty with a ``side`` disagreement when they differ. A number one reader prints stays in the
    table untrusted (``readers: 1``, ``confidence: low``; the engine ignores it) and becomes a
    ``number`` disagreement. Numbers are compared normalised, so ``569/C`` and ``569C`` agree.
    """
    keys = list(surveys) if surveys is not None else sorted(set(reading_a) | set(reading_b), key=paths.survey_sort_key)
    table: Dict[str, List[dict]] = {}
    disagreements: List[Disagreement] = []
    stats: List[SheetStats] = []
    for s in keys:
        a, b = _by_number(reading_a.get(s, [])), _by_number(reading_b.get(s, []))
        own = neighbours.normalise(s)
        numbers = sorted((set(a) | set(b)) - {own}, key=paths.survey_sort_key)
        rows: List[dict] = []
        checked = agreed = 0
        for n in numbers:
            ea, eb = a.get(n, []), b.get(n, [])
            if ea and eb:
                sides_a = {e.side for e in ea if e.side}
                sides_b = {e.side for e in eb if e.side}
                common = sorted(sides_a & sides_b)
                checked += 1
                if common:
                    agreed += 1
                    for side in common:
                        rows.append({"number": n, "side": side, "readers": 2, "confidence": "high", "source": SOURCE_READER})
                else:
                    rows.append({"number": n, "side": "", "readers": 2, "confidence": "high", "source": SOURCE_READER})
                    disagreements.append(Disagreement(village, s, "side", n, name_a, ea[0].number, "/".join(sorted(sides_a)),
                                                      name_b, eb[0].number, "/".join(sorted(sides_b)),
                                                      note="number agreed, side differs"))
            else:
                e, who = (ea[0], name_a) if ea else (eb[0], name_b)
                rows.append({"number": n, "side": e.side, "readers": 1, "confidence": "low", "source": SOURCE_READER})
                disagreements.append(Disagreement(village, s, "number", n,
                                                  name_a, ea[0].number if ea else "", ea[0].side if ea else "",
                                                  name_b, eb[0].number if eb else "", eb[0].side if eb else "",
                                                  note="read by %s only" % who))
        table[s] = rows
        stats.append(SheetStats(s, len(set(a) - {own}), len(set(b) - {own}), len(set(a) & set(b) - {own}),
                                len(numbers), checked, agreed))
    return table, disagreements, stats


def apply_resolutions(table: Dict[str, List[dict]], disagreements: List[Disagreement]) -> int:
    """Fold the analyst's ``resolution`` column back into the table; returns the count applied.

    ``number`` rows: yes keeps the number as verified (``readers: 2``), no drops it.
    ``side`` rows: a compass value becomes the side. Every field set here carries
    ``source = hand`` and ``resolved_by = hand``; reader accuracy is computed before this step
    and never counts these.
    """
    applied = 0
    for d in disagreements:
        res = (d.resolution or "").strip().upper()
        if not res:
            continue
        rows = table.get(d.survey, [])
        if d.field == "number":
            if res.lower() in RESOLVE_YES:
                for r in rows:
                    if r["number"] == d.number:
                        r["readers"], r["confidence"] = 2, "high"
                        r["source"], r["resolved_by"] = SOURCE_HAND, SOURCE_HAND
                d.resolved_by = SOURCE_HAND
                applied += 1
            elif res.lower() in RESOLVE_NO:
                table[d.survey] = [r for r in rows if r["number"] != d.number]
                d.resolved_by = SOURCE_HAND
                applied += 1
            else:
                log.warning("%s %s: resolution %r for a number row is not yes/no", d.survey, d.number, d.resolution)
        elif d.field == "side":
            if res in SIDES:
                for r in rows:
                    if r["number"] == d.number:
                        r["side"] = res
                        r["source"], r["resolved_by"] = SOURCE_HAND, SOURCE_HAND
                d.resolved_by = SOURCE_HAND
                applied += 1
            else:
                log.warning("%s %s: resolution %r is not a compass side", d.survey, d.number, d.resolution)
    return applied


# ----------------------------------------------------------------------------- files


def override_path(village: str) -> Path:
    return paths.vector_dir(village) / "neighbour_transcription" / ("nb_override_%s.json" % village)


def review_path(village: str) -> Path:
    return paths.vector_dir(village) / "neighbour_transcription" / ("review_%s.csv" % village)


def report_path(village: str) -> Path:
    return paths.vector_dir(village) / "neighbour_transcription" / ("agreement_%s.md" % village)


def seeds_path(village: str) -> Path:
    return paths.vector_dir(village) / "neighbour_transcription" / ("seeds_%s.json" % village)


def write_override(village: str, table: Dict[str, List[dict]]) -> Path:
    p = override_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(table, indent=1), encoding="utf-8")
    neighbours._CACHE.pop(village, None)
    return p


def read_review(path: Path) -> Dict[Tuple[str, str, str], Dict[str, str]]:
    """Existing review rows keyed by (survey, field, number), so resolutions survive a re-run."""
    if not Path(path).exists():
        return {}
    with Path(path).open(newline="", encoding="utf-8") as fh:
        return {(r["survey"], r["field"], r["number"]): r for r in csv.DictReader(fh)}


def carry_resolutions(disagreements: List[Disagreement], path: Path) -> int:
    old = read_review(path)
    n = 0
    for d in disagreements:
        prev = old.get(d.key())
        if prev and prev.get("resolution"):
            d.resolution = prev["resolution"]
            n += 1
    return n


def write_review_csv(village: str, disagreements: List[Disagreement]) -> Path:
    p = review_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=REVIEW_COLUMNS)
        w.writeheader()
        for d in disagreements:
            row = {k: getattr(d, k) for k in REVIEW_COLUMNS}
            row["resolved_by"] = SOURCE_HAND if (d.resolution or "").strip() else ""
            w.writerow(row)
    return p


def agreement_totals(stats: List[SheetStats]) -> Dict[str, object]:
    both = sum(s.n_both for s in stats)
    union = sum(s.n_union for s in stats)
    n_a, n_b = sum(s.n_a for s in stats), sum(s.n_b for s in stats)
    checked, agreed = sum(s.sides_checked for s in stats), sum(s.sides_agreed for s in stats)
    return {"sheets": len(stats), "labels_a": n_a, "labels_b": n_b, "numbers_both": both, "numbers_union": union,
            "number_jaccard": round(both / union, 3) if union else None,
            "number_dice": round(2 * both / (n_a + n_b), 3) if (n_a + n_b) else None,
            "sides_checked": checked, "sides_agreed": agreed,
            "side_agreement": round(agreed / checked, 3) if checked else None,
            "sheets_fully_agreed": sum(1 for s in stats if s.n_union and s.n_both == s.n_union and s.sides_agreed == s.sides_checked)}


def sample_agreed(stats: List[SheetStats], k: int, seed: int) -> List[str]:
    """``k`` random sheets on which the readers agreed completely, for the analyst's check."""
    full = [s.survey for s in stats if s.n_union and s.n_both == s.n_union and s.sides_agreed == s.sides_checked]
    rng = random.Random(seed)
    return sorted(rng.sample(full, min(k, len(full))), key=paths.survey_sort_key)


def write_agreement_report(village: str, stats: List[SheetStats], disagreements: List[Disagreement],
                           sample: List[str], readers: Sequence[str], extra: Optional[str] = None) -> Path:
    t = agreement_totals(stats)
    lines = ["# Two-reader transcription: %s" % village, "",
             datetime.datetime.now().isoformat(timespec="seconds"), "",
             "Readers: %s. Both are the same vision model in separate contexts, so agreement between "
             "them overstates accuracy against the PDF; the analyst's sample check below is the accuracy "
             "measure. Trusted entries are those both readers printed (readers = 2)." % ", ".join(readers), "",
             "## Agreement", "",
             "| measure | value |", "|---|---|"]
    for k, val in t.items():
        lines.append("| %s | %s |" % (k, val))
    lines += ["", "## Disagreements for review (%d)" % len(disagreements), "",
              "Listed in the review CSV with a `resolution` column: for a `number` row write yes or no, "
              "for a `side` row write the compass side.", "",
              "| sheet | field | number | reader A | reader B | note |", "|---|---|---|---|---|---|"]
    for d in disagreements:
        lines.append("| %s | %s | %s | %s %s | %s %s | %s |" % (d.survey, d.field, d.number, d.printed_a, d.side_a, d.printed_b, d.side_b, d.note))
    lines += ["", "## Sample of fully agreed sheets to check against the PDF (%d)" % len(sample), "",
              ", ".join(sample) or "(none)", "",
              "## Per sheet", "", "| sheet | labels A | labels B | both | union | jaccard | sides agreed / checked |", "|---|---|---|---|---|---|---|"]
    for s in stats:
        lines.append("| %s | %d | %d | %d | %d | %s | %d / %d |" % (s.survey, s.n_a, s.n_b, s.n_both, s.n_union, s.jaccard, s.sides_agreed, s.sides_checked))
    if extra:
        lines += ["", extra]
    p = report_path(village)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


# ----------------------------------------------------------------------------- graph and seeds


def neighbour_graph(table: Dict[str, List[dict]], surveys: Iterable[str]) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    """Undirected graph over ``surveys`` from trusted entries: an edge when either sheet prints
    the other. Second result: printed numbers that are not in ``surveys`` (no sheet), per sheet."""
    keep = {neighbours.normalise(s) or str(s).upper(): str(s) for s in surveys}
    graph: Dict[str, Set[str]] = {s: set() for s in keep.values()}
    outside: Dict[str, Set[str]] = {}
    for s, rows in table.items():
        s_key = neighbours.normalise(s) or str(s).upper()
        if s_key not in keep:
            continue
        me = keep[s_key]
        for r in rows:
            if int(r.get("readers", 1)) < 2 and r.get("confidence") != "high":
                continue
            n = r["number"]
            if n in keep:
                other = keep[n]
                if other != me:
                    graph[me].add(other)
                    graph[other].add(me)
            else:
                outside.setdefault(me, set()).add(n)
    return graph, outside


def corrected_sheets(village: str) -> Set[str]:
    """Surveys whose OUTLINE was corrected by hand (sheet_corrections.csv rows other than plot
    re-conversions). They are never seeds and never labels."""
    p = paths.vector_dir(village) / "sheet_corrections.csv"
    if not p.exists():
        return set()
    with p.open(newline="", encoding="utf-8") as fh:
        return {str(r.get("survey")) for r in csv.DictReader(fh) if (r.get("edge") or "") != "plots"}


def dominating_set(graph: Dict[str, Set[str]], fixed: FrozenSet[str] = frozenset(),
                   prefer: FrozenSet[str] = frozenset(), exclude: FrozenSet[str] = frozenset()) -> List[str]:
    """Greedy dominating set: every node ends up a seed or adjacent to one.

    ``fixed`` nodes are seeds already (nothing to choose, they cost nothing); ``prefer`` breaks
    ties towards nodes the analyst has already placed; ``exclude`` nodes may never be seeds
    (hand-corrected sheets) but still need covering. Deterministic: ties then fall to the
    survey sort order.
    """
    nodes = set(graph)
    covered: Set[str] = set()
    for f in fixed:
        if f in graph and f not in exclude:
            covered |= {f} | graph[f]
    chosen: List[str] = []
    while nodes - covered:
        uncovered = nodes - covered
        candidates = nodes - set(chosen) - set(fixed) - set(exclude)
        if not candidates:
            break
        # most new coverage first, then a node already placed, then survey order
        best = min(candidates, key=lambda n: (-len(({n} | graph[n]) & uncovered), n not in prefer, paths.survey_sort_key(n)))
        gain = len(({best} | graph[best]) & uncovered)
        if gain == 0:
            # isolated nodes: each is its own seed
            for n in sorted(uncovered, key=paths.survey_sort_key):
                chosen.append(n)
                covered.add(n)
            break
        chosen.append(best)
        covered |= {best} | graph[best]
    return chosen


def seed_plan(village: str, table: Dict[str, List[dict]], surveys: Iterable[str], placed: Iterable[str]) -> Dict[str, object]:
    """Both readings of the seed question: from scratch (prefer his placements) and additional to
    what he has placed. Seeds are anchors and never scored; labelled parcels are his placements
    that are not seeds."""
    surveys = list(surveys)
    placed = frozenset(str(p) for p in placed)
    excluded = frozenset(corrected_sheets(village))
    graph, outside = neighbour_graph(table, surveys)
    scratch = dominating_set(graph, prefer=placed - excluded, exclude=excluded)
    additional = dominating_set(graph, fixed=placed - excluded, exclude=excluded)
    isolated = sorted((n for n, nb in graph.items() if not nb), key=paths.survey_sort_key)
    return {
        "village": village, "sheets": len(surveys), "edges": sum(len(v) for v in graph.values()) // 2,
        "isolated_sheets": isolated,
        "placed": sorted(placed, key=paths.survey_sort_key),
        "excluded_from_seeds_and_labels": sorted(excluded, key=paths.survey_sort_key),
        "from_scratch": {"seeds": sorted(scratch, key=paths.survey_sort_key), "n": len(scratch),
                         "already_placed": sorted(set(scratch) & placed, key=paths.survey_sort_key),
                         "to_place": sorted(set(scratch) - placed, key=paths.survey_sort_key),
                         "labelled_parcels": sorted(placed - set(scratch) - excluded, key=paths.survey_sort_key)},
        "additional": {"seeds": sorted(additional, key=paths.survey_sort_key), "n": len(additional),
                       "labelled_parcels": sorted(placed - excluded, key=paths.survey_sort_key)},
        "printed_without_sheet": {k: sorted(v, key=paths.survey_sort_key) for k, v in sorted(outside.items(), key=lambda kv: paths.survey_sort_key(kv[0]))},
    }


# ----------------------------------------------------------------------------- orchestration


def run(village: str, cfg: Optional[configmod.Config] = None, placed: Optional[Iterable[str]] = None) -> Dict[str, object]:
    """Merge the two readings of a village and write everything the analyst and the engine need."""
    cfg = cfg or configmod.load()
    specs = list(cfg.transcription.readers)
    readings = load_readings(village, cfg)
    ra, rb = readings[specs[0]], readings[specs[1]]
    surveys = paths.surveys_with_sheets(village)
    missing = [s for s in surveys if s not in ra or s not in rb]
    table, disagreements, stats = merge(ra, rb, surveys, village, specs[0], specs[1])
    carried = carry_resolutions(disagreements, review_path(village))
    applied = apply_resolutions(table, disagreements)
    p_override = write_override(village, table)
    p_review = write_review_csv(village, disagreements)
    sample = sample_agreed(stats, cfg.transcription.sample_size, cfg.transcription.sample_seed)
    if placed is None:
        from . import anchors
        placed = list(anchors.load_anchors(village, set()))
    plan = seed_plan(village, table, surveys, placed)
    seeds_path(village).write_text(json.dumps(plan, indent=1), encoding="utf-8")
    n_hand = sum(1 for d in disagreements if (d.resolution or "").strip())
    extra = ("## Sources\n\nEvery table entry carries `source`: `reader` (both readers agreed) or `hand` (set by the "
             "analyst through the review CSV, `resolved_by = hand`). The agreement figures above are computed before "
             "resolutions and therefore exclude hand-resolved rows; %d row(s) are hand-resolved so far.\n\n"
             "## Seeds\n\nFrom scratch: %d seeds (%d already placed, %d to place: %s); %d labelled parcels.\n"
             "Additional to the %d placed: %d more seeds: %s.\nIsolated sheets (no trusted neighbour with a sheet): %s.\n"
             "Excluded from seeds and labels (outline corrected by hand): %s.\n"
             % (n_hand, plan["from_scratch"]["n"], len(plan["from_scratch"]["already_placed"]), len(plan["from_scratch"]["to_place"]),
                ", ".join(plan["from_scratch"]["to_place"]) or "-", len(plan["from_scratch"]["labelled_parcels"]),
                len(plan["placed"]), plan["additional"]["n"], ", ".join(plan["additional"]["seeds"]) or "-",
                ", ".join(plan["isolated_sheets"]) or "-", ", ".join(plan["excluded_from_seeds_and_labels"]) or "-"))
    if missing:
        extra += "\nSheets missing from a reading: %s.\n" % ", ".join(missing)
    p_report = write_agreement_report(village, stats, disagreements, sample, specs, extra)
    totals = agreement_totals(stats)
    log.info("%s: %d sheets, number agreement %s, side agreement %s, %d disagreements, %d resolutions applied",
             village, len(stats), totals["number_jaccard"], totals["side_agreement"], len(disagreements), applied)
    return {"village": village, "sheets": len(stats), "missing": missing, "totals": totals,
            "disagreements": len(disagreements), "resolutions_carried": carried, "resolutions_applied": applied,
            "sample": sample, "override": str(p_override), "review": str(p_review), "report": str(p_report),
            "seeds": str(seeds_path(village)), "plan": plan}
