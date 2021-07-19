# Zig Minesweeper Solver Pieces

Zig code that may be used for solving minesweeper boards (possibly in combination with the rest of this package!).

Uses Zig 0.8.


## How to Use

- Fetch the git submodule dependencies (`git submodule update`).
- Build with `zig build` (from the `zig/` subdirectory)
- Run with `./zig-out/bin/zig-main`

For usage see `zig-main --help`.

Example usage:
```
$./zig-out/bin/zig-main 3 <<EOF
>
> #  #  1
> 1  3  *
> 0  2  *
>
>EOF
# # 1
1 3 *
. 2 *
```

This shows the program reading in from stdin (use `-f` to read from a file), parsing the input, and outputting the board in canonical form.

The following cell representations are used:
- `#` for an unclicked cell
- `<N>` where `N=0,1,2,...` is a number shown in a cell
- `.` as an alternative to `0` (since the number 0 is not normally shown)
- `*` to represent a single mine (may be a revealed mine or a flag)
- `*<N>` where `N=1,2,...` is the number of mines
