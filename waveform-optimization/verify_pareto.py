"""Verify the Pareto and Knee Point calculations from the waveform optimization CSV files."""

import csv
from collections import Counter
from decimal import Decimal
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent


def read_rows(filename):
    with (DATA_DIR / filename).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def objective_vector(row):
    return (
        -Decimal(row["J_relative_to_DS_end"]),
        Decimal(row["C_relative_to_DS_end"]),
        Decimal(row["R_tail_percent"]),
    )


def dominates(first, second):
    return all(a <= b for a, b in zip(first, second)) and any(
        a < b for a, b in zip(first, second)
    )


def knee_distances(objectives):
    columns = list(zip(*objectives.values()))
    lower = [min(column) for column in columns]
    span = [max(column) - minimum for column, minimum in zip(columns, lower)]
    return {
        candidate: sum(
            ((value - minimum) / width) ** 2
            for value, minimum, width in zip(vector, lower, span)
        ).sqrt()
        for candidate, vector in objectives.items()
    }


def main():
    rows = read_rows("candidate_objectives.csv")
    shapes = {row["candidate_id"]: row for row in rows}
    if len(shapes) != len(rows):
        raise ValueError("The objective table contains duplicate candidate IDs.")
    objectives = {candidate: objective_vector(row) for candidate, row in shapes.items()}
    front = [
        candidate for candidate, vector in objectives.items()
        if not any(dominates(other, vector) for other in objectives.values())
    ]
    distances = knee_distances(objectives)

    for candidate, row in shapes.items():
        flag = "YES" if candidate in front else "NO"
        if row["pareto_status"] != flag:
            raise ValueError(f"Pareto status differs for {candidate}.")
        example = row["example_dominator_id"]
        if example and not dominates(objectives[example], objectives[candidate]):
            raise ValueError(f"{example} does not dominate {candidate}.")
        saved_distance = Decimal(row["knee_distance_full_shape_set"])
        if abs(distances[candidate] - saved_distance) > Decimal("1e-8"):
            raise ValueError(f"Knee distance differs for {candidate} at CSV precision.")

    mapping = read_rows("parameter_shape_mapping.csv")
    if len({row["raw_id"] for row in mapping}) != len(mapping):
        raise ValueError("The mapping table contains duplicate raw IDs.")
    if {row["shape_id"] for row in mapping} != set(shapes):
        raise ValueError("The mapping and objective tables have different shape IDs.")
    for row in mapping:
        shape = shapes[row["shape_id"]]
        parameter = shape["shape_parameter"]
        if row["family"] != shape["family"] or Decimal(row[parameter]) != Decimal(shape["shape_value"]):
            raise ValueError(f"Shape parameters differ for {row['raw_id']}.")

    knee = min(front, key=distances.__getitem__)
    families = Counter(shapes[candidate]["family"] for candidate in front)
    print(f"{len(mapping):,} nominal parameter combinations map to {len(shapes)} shapes.")
    print(f"Nondominated shapes: {len(front)} ({', '.join(f'{name}: {n}' for name, n in sorted(families.items()))})")
    print(f"Nominal Knee Point: {knee} (normalized distance {distances[knee]:.8f})")
    for figure in (5, 6):
        count = sum(row[f"shown_in_Fig{figure}"] == "YES" for row in rows)
        print(f"Figure {figure}: {count} plotted shape records")


if __name__ == "__main__":
    main()
