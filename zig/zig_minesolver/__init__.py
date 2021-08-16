__all__ = ("get_board_probs",)

import pathlib
import subprocess
from typing import List, Optional


THIS_DIR = pathlib.Path(__file__).parent


def get_board_probs(
    board: str,
    *,
    mines: Optional[int] = None,
    density: Optional[float] = None,
    per_cell: int = 1,
) -> List[List[float]]:
    """
    Get probabilities for a minesweeper board.

    Takes a board in string format, where:
     - `#` represents an unclicked cell
     - `<N>` where `N=0,1,2,...` is a number shown in a cell
     - `.` as an alternative to `0` (since the number 0 is not normally shown)
     - `*` to represent a single mine (may be a revealed mine or a flag)
     - `*<N>` where `N=1,2,...` is the number of mines

    Returns a 2-D list representation of the board, where each cell is the
    probability of the cell being a mine, between 0 and 1.

    Either a number of mines or a density of mines must be given.
    """
    cmd = [str(THIS_DIR / "zig-main"), "--per-cell", str(per_cell)]
    if mines is not None:
        cmd += ["--mines", str(mines)]
    if density is not None:
        cmd += ["--infinite-density", str(density)]
    proc = subprocess.run(
        cmd,
        input=board,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=5,
        check=True,
    )
    lines = proc.stdout.splitlines()
    start_idx = lines.index("Probabilities:") + 1
    return [[float(x) for x in L.split()] for L in lines[start_idx:]]
