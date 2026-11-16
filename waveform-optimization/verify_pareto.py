"""Check the Pareto front and nominal Knee Point for the candidate table."""

import csv
from collections import Counter
from decimal import Decimal
from pathlib import Path


DATA_DIR = Path(__file__).resolve().parent


def read_rows(filename):
    with (DATA_DIR / filename).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def objective_vector(row):
    # Negating J puts all three objectives in minimization form.
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
    # The min-max bounds come from the full table, including dominated points.
    bounds = [(min(column), max(column)) for column in zip(*objectives.values())]
    distances = {}
    for candidate_id, vector in objectives.items():
        squared_distance = sum(
            ((value - low) / (high - low)) ** 2
            for value, (low, high) in zip(vector, bounds)
        )
        distances[candidate_id] = squared_distance.sqrt()
    return distances


def main():
    rows = read_rows("candidate_objectives.csv")
    candidates = {row["candidate_id"]: row for row in rows}
    if len(candidates) != len(rows):
        raise ValueError("Duplicate candidate IDs in candidate_objectives.csv.")

    objectives = {candidate_id: objective_vector(row) for candidate_id, row in candidates.items()}
    front = [
        candidate_id for candidate_id, vector in objectives.items()
        if not any(dominates(other, vector) for other in objectives.values())
    ]
    distances = knee_distances(objectives)

    distance_tolerance = Decimal("1e-8")
    for candidate_id, row in candidates.items():
        status = "YES" if candidate_id in front else "NO"
        if row["pareto_status"] != status:
            raise ValueError(f"Expected Pareto status {status} for {candidate_id}.")
        dominator_id = row["example_dominator_id"]
        if dominator_id and not dominates(objectives[dominator_id], objectives[candidate_id]):
            raise ValueError(f"{dominator_id} does not dominate {candidate_id}.")
        saved_distance = Decimal(row["knee_distance_full_shape_set"])
        if abs(distances[candidate_id] - saved_distance) > distance_tolerance:
            raise ValueError(f"Knee distance differs for {candidate_id} at CSV precision.")

    mapping = read_rows("parameter_shape_mapping.csv")
    if len({row["raw_id"] for row in mapping}) != len(mapping):
        raise ValueError("Duplicate raw IDs in parameter_shape_mapping.csv.")
    if {row["shape_id"] for row in mapping} != set(candidates):
        raise ValueError("The shape IDs in the mapping must match the candidate IDs in the objective table.")
    for row in mapping:
        candidate = candidates[row["shape_id"]]
        parameter = candidate["shape_parameter"]
        if (
            row["family"] != candidate["family"]
            or Decimal(row[parameter]) != Decimal(candidate["shape_value"])
        ):
            raise ValueError(f"Shape parameters differ for {row['raw_id']}.")

    knee = min(front, key=distances.__getitem__)
    families = Counter(candidates[candidate_id]["family"] for candidate_id in front)
    counts = ", ".join(f"{name}: {count}" for name, count in sorted(families.items()))
    print(f"{len(mapping):,} nominal parameter combinations map to {len(candidates)} candidate records.")
    print(f"Nondominated shapes: {len(front)} ({counts})")
    print(f"Nominal Knee Point: {knee} (normalized distance {distances[knee]:.8f})")
    for figure in (5, 6):
        count = sum(row[f"shown_in_Fig{figure}"] == "YES" for row in rows)
        print(f"Figure {figure}: {count} plotted shape records")


if __name__ == "__main__":
    main()
