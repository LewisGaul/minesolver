import subprocess
from typing import List


def get_board_probs(board: str, mines: int, *, per_cell: int = 1) -> List[List[float]]:
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
    """
    proc = subprocess.run(
        ["zig-main", str(mines), "--per-cell", str(per_cell)],
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
