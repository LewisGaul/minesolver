# Zig Minesweeper Solver

Zig code that may be used for solving minesweeper boards (possibly in combination with the rest of this package!).

Uses Zig 0.11.


## How to Use

- Fetch the git submodule dependencies (`git submodule update`).
- Build with `zig build` (from the `zig/` subdirectory)
- Run with `./zig-out/bin/zig-main`

For usage see `zig-main --help`.

The following cell representations are used for input boards:
- `#` for an unclicked cell
- `<N>` where `N=0,1,2,...` is a number shown in a cell
- `.` as an alternative to `0` (since the number 0 is not normally shown)
- `*` to represent a single mine (may be a revealed mine or a flag)
- `*<N>` where `N=1,2,...` is the number of mines


Example usages:
```
$time ./zig-out/bin/zig-main -f example2.txt --mines 8 2>/dev/null

Board:
# 2 # # #
# # # # #
# 3 # # #
# 2 # 4 #
# # # # #

Matrix:
1 1 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 | 2
0 0 0 0 1 1 1 0 0 1 1 0 0 1 1 0 0 0 0 0 0 | 3
0 0 0 0 0 0 0 0 0 1 1 0 0 1 1 0 1 1 1 0 0 | 2
0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 1 0 0 1 1 1 | 4
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 | 8

RREF matrix:
 1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0  0 |  1
 0  0  1  1  0  0  0  1  1  0  0  1  1  0  0  1  0  0  0  1  1 |  4
 0  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0 -1 -1 -1  0  0 |  1
 0  0  0  0  0  0  0  0  0  1  0 -1 -1  1  0 -1  1  1  0 -1 -1 | -2
 0  0  0  0  0  0  0  0  0  0  1  1  1  0  1  1  0  0  1  1  1 |  4

Solver matrix:
 1  0  0  0  0  0  1  1 |  1
 0  1  0  0  0  1  0  0 |  4
 0  0  1  0  0  0 -1 -1 |  1
 0  0  0  1  0 -1  1  0 | -2
 0  0  0  0  1  1  0  1 |  4

Solver groups:
0: { 0, 1 }
1: { 2, 3, 7, 8 }
2: { 4, 5, 6 }
3: { 9, 13 }
4: { 10, 14 }
5: { 11, 12, 15, 19, 20 }
6: { 16, 17 }
7: { 18 }

Mine configurations:
0: { 1, 2, 1, 0, 2, 2, 0, 0 }
1: { 1, 1, 1, 1, 1, 3, 0, 0 }
2: { 1, 0, 1, 2, 0, 4, 0, 0 }
3: { 0, 1, 2, 0, 1, 3, 1, 0 }
4: { 0, 0, 2, 1, 0, 4, 1, 0 }
5: { 0, 2, 2, 0, 1, 2, 0, 1 }
6: { 0, 1, 2, 1, 0, 3, 0, 1 }

Probabilities:
0.2711 0.0000 0.2711 0.3133 0.3133
0.4859 0.4859 0.4859 0.3133 0.3133
0.2651 0.0000 0.5060 0.5494 0.5494
0.2651 0.0000 0.5060 0.0000 0.5494
0.1084 0.1084 0.2410 0.5494 0.5494

real    0m0.007s
user    0m0.000s
sys     0m0.000s
```


```
$time ./zig-out/bin/zig-main -f example3.txt --mines 99 --per-cell 3 >/dev/null
 0.005 [ WARN]: Omitting large full matrix output
 0.006 [ WARN]: Omitting large RREF matrix output
 0.006 [ INFO]: Initialising solver with 30 x 16 board
 0.009 [ INFO]: Initial matrix is 388 x 75
 0.017 [ INFO]: Reduced matrix to groups, 52 columns
 0.017 [ INFO]: Reduced matrix to RREF
 0.017 [ INFO]: Removed zero-rows from matrix, 42 rows
 0.029 [ INFO]: Categorised columns: 42 fixed, 9 free
 0.029 [ INFO]: Using free val maximums: { 3, 3, 2, 3, 1, 3, 1, 3, 0 } (12288 combinations)
 0.062 [ INFO]: Found 60 mine configurations
 0.104 [ INFO]: Finished

real    0m0.117s
user    0m0.047s
sys     0m0.063s
```
