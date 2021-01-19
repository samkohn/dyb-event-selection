import contextlib
import sqlite3

__all__ = [
        "ALL_DETECTORS",
        "AD_DETECTORS",
        "WP_DETECTORS",
        "dets_for",
        ]
ALL_DETECTORS = [0, 1, 2, 3, 4, 5, 6, 7, 8]
AD_DETECTORS = [1, 2, 3, 4]
WP_DETECTORS = [5, 6]

near_ads = [(1, 1), (1, 2), (2, 1), (2, 2)]
far_ads = [(3, 1), (3, 2), (3, 3), (3, 4)]
all_ads = near_ads + far_ads

def phase_for_run(runno):
    if 21221 <= runno <= 26694:
        return 1
    if 34523 <= runno <= 67012:
        return 2
    if runno >= 67625:
        return 3
    raise ValueError(f"Nonsensical run number: {runno}")

def phase_for_day(day):
    if 0 <= day <= 217:
        return 1
    if 300 <= day <= 1824:
        return 2
    if 1860 <= day <= 2076:
        return 3
    raise ValueError(f'Day {day} is not part of P17B')

def dets_for_phase(site, phase):
    if site == 1:
        return [2] if phase == 3 else [1, 2]
    if site == 2:
        return [1] if phase == 1 else [1, 2]
    if site == 3:
        return [1, 2, 3] if phase == 1 else [1, 2, 3, 4]
    raise ValueError("Invalid site")

def dets_for(site, runno=None):
    if runno is None:
        phase = 2  # 8-ad, i.e. return all ADs for site
    else:
        phase = phase_for_run(runno)
    return dets_for_phase(site, phase)


@contextlib.contextmanager
def get_db(db, *args, **kwargs):
    """Return a context manager that closes the db automatically.

    This is stupid because the "with sqlite3.connect/Connection"
    will not close the db connection, and "with closing(sqlite3.connect)"
    will not commit the transaction. So you need both.
    """
    with contextlib.closing(sqlite3.connect(db, *args, **kwargs)) as conn:
        with conn:
            yield conn
