"""Every test must pass in at least one CI job.

Each job skips the tests for tools it doesn't have (the CAD job has no DOLFINx,
the FEniCSx job has no CadQuery or ccx). This reads both jobs' JUnit reports and
fails if any test was skipped or failed everywhere.

    python scripts/check_no_skips.py out/junit-cad.xml out/junit-fenicsx.xml
"""
import sys
import xml.etree.ElementTree as ET


def outcomes(path: str) -> dict[str, str]:
    res = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        if case.find("skipped") is not None:
            res[key] = "skipped"
        elif case.find("failure") is not None or case.find("error") is not None:
            res[key] = "failed"
        else:
            res[key] = "passed"
    return res


def main(paths: list[str]) -> int:
    runs = [outcomes(p) for p in paths]
    tests = sorted(set().union(*runs))
    passed = [t for t in tests if any(r.get(t) == "passed" for r in runs)]
    missing = [t for t in tests if t not in passed]
    for p, r in zip(paths, runs):
        n = {k: sum(1 for v in r.values() if v == k) for k in ("passed", "skipped", "failed")}
        print(f"{p}: {n['passed']} passed, {n['skipped']} skipped, {n['failed']} failed")
    print(f"{len(tests)} tests in all; {len(passed)} passed in at least one job")
    for t in missing:
        print(f"NOT RUN ANYWHERE: {t}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
