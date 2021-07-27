# Solver Design

This document contains an implementation-level design of this minesweeper solver.


## Probability Calculation

### Matrix approach

The first step in the process is to transform the board into a structure that can be worked with more easily. This makes sense to be done by transforming to matrix form (since a minesweeper board is effectively just a set of simultaneous equations). An alternative approach is sketched out in a section below.

- Each row of the matrix corresponds to one of the simultaneous equations, which come from each visible number on the board.
- There is an additional final row corresponding to the equation for the total number of mines in the board.
- Each column corresponds to an unclicked cell, where the value in the matrix is '1' if the unclicked cell is a neighbour of the row's displayed number, and '0' otherwise.
- There is a single column on the right corresponding to the RHS of the simultaneous equations, i.e. the value of the number shown in the cell corresponding to that row.

To find the possible configurations of mines in the board, we just need to find all non-negative integer solutions (also bounded by the maximum number of mines per cell).

The matrix can be reduced to RREF (reduced-row echelon form). This is not strictly necessary, but can help with reading off solutions, and can serve as an optimisation for finding solutions. Since the initial matrix contains only integers, the reduction can be performed such that the resulting matrix also contains only integers (positive or negative) - this is at the cost of being able to guarantee each leading row value is a '1', but this is inconsequential.

Example:
```
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
```


### Finding equivalence groups

Note that this part of the process is not strictly necessary, but does provide some insight into the board structure, and can serve as an optimisation. It does involve a bit of careful bookkeeping though.

When using the matrix approach this is as simple as finding identical columns - these correspond to unclicked cells that are in the same equivalence group, and they can be collapsed into a single column. In this case each column of the resulting matrix corresponds to a group (which in turn corresponds to a number of cells).

Example:

```
Board:
# 2 # # #
# # # # #
# 3 # # #
# 2 # 4 #
# # # # #

Matrix:
1 0 1 0 0 0 0 0 | 2
0 0 1 1 1 0 0 0 | 3
0 0 0 1 1 0 1 1 | 2
0 0 0 0 1 1 0 1 | 4
1 1 1 1 1 1 1 1 | 8

RREF matrix:
 1  0  0  0  0  0  1  1 |  1
 0  1  0  0  0  1  0  0 |  4
 0  0  1  0  0  0 -1 -1 |  1
 0  0  0  1  0 -1  1  0 | -2
 0  0  0  0  1  1  0  1 |  4

Groups (<group number>: <unclicked cell numbers>):
1: 1, 2
2: 3, 4, 8, 9   <- This is the 'outer' group, not next to any numbers
3: 5, 6, 7
4: 10, 14
5: 11, 16
6: 12, 13, 15, 20, 21
7: 17, 18
8: 19
```


### Finding mine configurations

As mentioned above, this just comes down to finding all non-negative integer solutions to the matrix equation (also bounded by maximum mines per cell).

For matrix equations over the real numbers there are three possibilities, which can immediately be determined when the matrix is in RREF.
- One unique solution - the top of the LHS of the RREF matrix is a square identity matrix and all rows below are zero (only need to check the first row below on the RHS in RREF) - the solution can be trivially read off.
- No solutions - there is a row with all zeros on the LHS and a non-zero value on the RHS.
- Infinite solutions - there are no rows with all zeros on the LHS and a non-zero value on the RHS, and there are less non-zero rows than there are columns.

These rules are still useful for determining constrained integer solutions, although extra consideration is needed. Turning around the cases above:
- After removing any zero-rows, the LHS of the RREF matrix is a square identity matrix - the unique potential solution can be read off and checked against the bound constraints (e.g. all values must be non-negative).
- There is at least one row with all zeros on the LHS and a non-zero value on the RHS - invalid board (no solutions).
- There are no rows with all zeros on the LHS and a non-zero value on the RHS, and there are less non-zero rows than there are columns - further work is needed to find valid solutions.

In most cases we expect a minesweeper board to fall into the latter case. We discuss the process for finding solutions in this case below.

Any zero-rows can be removed - these correspond to there being multiple displayed board numbers that provide exactly the same information, i.e. redundant information.

Columns of the RREF matrix can be classified into two cases:
- A column that contains a leading-edge entry (i.e. the first non-zero value in a row, with all other values in the column being zero) - we say these columns correspond to 'fixed' variables.
- Otherwise (if there are multiple non-zero values in the column, or if the non-zero value is not a leading-edge entry) the column corresponds to a 'free' variable.

The purpose of this is that you can write solutions for all fixed variables in terms of the free variables. For example, using the RREF matrix above, with the columns corresponding to groups `g1`, `g2`, ...:
```
 1  0  0  0  0  0  1  1 |  1
 0  1  0  0  0  1  0  0 |  4
 0  0  1  0  0  0 -1 -1 |  1
 0  0  0  1  0 -1  1  0 | -2
 0  0  0  0  1  1  0  1 |  4

Fixed variables: g1, g2, g3, g4, g5
Free variables: g6, g7, g8

g1 =  1 - g7 - g8
g2 =  4 - g6
g3 =  1 + g7 + g8
g4 = -2 + g6 - g7
g5 =  4 - g6 - g8
```

The problem then becomes a matter of finding all valid values for the free variables, with the constraint being: `0 <= g<n> <= |G<n>| * <max mines per cell>`, where `g<n>` is the number of mines in the n'th group and `|G<n>|` is the size of the n'th group.

A basic approach is to start with the free variables at 0, and increase them until any of the fixed variable values become invalid (or the max for the free variable is reached). Note that in practice it is more common for the non-negative constraint to be stronger when the maximum number of mines per cell is greater than 1. If we use the fact that the fixed variables are non-negative, we can write down an inequality from the above equations:
```
| g1 |   |  1 |   |  0  1  1 |             | 0 |
| g2 |   |  4 |   |  1  0  0 |   | g6 |    | 0 |
| g3 | = |  1 | - |  0 -1 -1 | . | g7 | >= | 0 |
| g4 |   | -2 |   | -1  1  0 |   | g8 |    | 0 |
| g5 |   |  4 |   |  1  0  1 |             | 0 |

|  0  1  1 |             |  1 |
|  1  0  0 |   | g6 |    |  4 |
|  0 -1 -1 | . | g7 | <= |  1 |
| -1  1  0 |   | g8 |    | -2 |
|  1  0  1 |             |  4 |
```

This can be used to check whether values chosen for the free variables are potential solutions, but we can also use it to set some tigher upper bounds on the free variables, using only the non-negativity property of all of the variables.

All of the rows that contain no negative values give bounds when you consider that the free variables are non-negative:
```
Row 1: g7 + g8 <= 1  =>  g7 <= 1, g8 <= 1
Row 2: g6 <= 4
Row 5: g6 + g8 <= 4  =>  g6 <= 4, g8 <= 4
```

At this point we have the following intervals to consider (assuming max per cell is not a limiting factor - this can be applied as an additional constraint without much difficulty if needed):
```
0 <= g6 <= 4
0 <= g7 <= 1
0 <= g8 <= 1
```

This gives a total of `5*2*2=20` possibilities to check, whereas without narrowing down in this way, even if the max mines per cell was 1, there would be `6*3*2=36` possibilities to check.

It turns out that there are 7 valid mine configurations in this example:
```
(1, 2, 1, 0, 2, 2, 0, 0)
(1, 1, 1, 1, 1, 3, 0, 0)
(1, 0, 1, 2, 0, 4, 0, 0)
(0, 1, 2, 0, 1, 3, 1, 0)
(0, 0, 2, 1, 0, 4, 1, 0)
(0, 2, 2, 0, 1, 2, 0, 1)
(0, 1, 2, 1, 0, 3, 0, 1)
```


### Calculating the probabilities

TODO



## Sketch Implementation

A sketch of the implementation is given below, using a Rust-like syntax to give a view of relevant APIs and their usage.


### Grid and Board types

There's an underlying grid type to model the 2-dimensional grid. There are then two types wrapping this: the board that we're working out probabilities for and the grid of probabilities for the board. At the highest-level view this is enough to encapsulate what we want to achieve: to map an in-progress game board to a grid of probabilities.

```rust
struct Grid<T> {
    fn cell(coord: Coord) -> T;
    fn coords() -> Vec<Coord>;
    fn nbr_coords(coord: Coord) -> Vec<Coord>;
}

struct Board : Grid<BoardCellContents> {
    mines: u32,    // number of mines
    per_cell: u32, // max mines per cell
    // 'Inherit' Grid methods
}

type ProbabilityGrid = Grid<f32>;
```


### Board reduction

TODO


### Finding mine configurations

Now the most complex (and compuatationally intensive) step - we need to take this reduced representation and find all possible ways to place mines in the equivalence groups of unclicked cells.

TODO


### Calculating probabilities

Finally we need to convert these combinations into a grid of probabilities.

```rust
fn find_probabilities(combs: Vec<Vec<u32>>, groups: Vec<BoardGroup>, per_cell: u32) -> ProbabilityGrid {
    // TODO
}
```


## Alternative Approach

This approach avoids the use of matrices - in some ways it's a bit easier to relate the process back to the board structure, but in some ways it's also harder to keep track of things.

[Note: We're actually skipping an optional step here of reducing the board by applying solver logic, e.g. if a number '1' is only next to one unclicked cell.]

```rust
struct DecomposedBoard {
    mines: u32,
    // Note the per_cell restriction is encapsulated in the groups.
    groups: Vec<BoardGroup>,
    numbers: Vec<BoardNumber>,
    
    fn get_group(coord: Coord) -> Option<&BoardGroup>;
    fn get_number(coord: Coord) -> Option<&BoardNumber>;
}

struct BoardGroup {
    // Coordinates of the unclicked cells in the equivalence group.
    coords: Vec<Coord>,
    // Neighbouring cells displaying a number.
    nbr_numbers: Vec<&Number>,
    // Max number of mines that can be contained in the group.
    max_mines: u32,
}

struct BoardNumber {
    coord: Coord,            // coord of the number cell
    value: u32,              // value of the number
    nbr_groups: Vec<&Group>, // neighbouring groups
}
```

```rust
fn decompose_board(board: &Board) -> DecomposedBoard {
    let mut decomp = DecomposedBoard::new(board.mines);
    
    // This will be a mapping of unclicked cell coords to coords of
    // neighbouring numbers: {(x, y): [(a, b), (c, d), ...], ... }
    let mut unclicked_cell_nbr_nums = HashMap::new();
    
    for (num_coord, num_val) in board.iter_num_cells() {
        // We only care about numbers next to unclicked cells.
        for unclicked in board.unclicked_nbr_coords() {
            unclicked_cell_nbr_nums[unclicked].push(num_coord)
        }
    }

    // Remap to {[(a, b), ...]: [(x, y), ...], ... }
    let nums_to_groups = unclicked_cell_nbr_nums.group_by_value();
    for (num_coords, group_coords) in nums_to_groups {
        let mut group = BoardGroup::new(group_coords);
        for num_coord in num_coords {
            // Create BoardNumber if not yet created.
            let mut number: &BoardNumber = decomp.get_or_add_number(num_coord);
            // Store references between numbers and groups.
            group.numbers.push(number);
            number.groups.push(&group);
        }
        decomp.groups.push(group);
    }
    
    decomp
}
```


```rust
fn find_combinations(board: DecomposedBoard) -> Combinations {
    /*
     * Example of what we have here:
     *
     *   # 2 # # #       . 2 . . M .
     *   # # # # #       . M M . . .
     *   # 3 # # #       . 3 . M M .
     *   # 2 # 4 #       M 2 . 4 M .
     *   # # # # #       . . M . . .
     *
     * DecomposedBoard {
     *   mines: 8,
     *   numbers: [
     *     {(1,0), 2, groups: [0, 1]},
     *     {(1,2), 3, groups: [1, 2, 4]},
     *     {(1,3), 2, groups: [2, 3, 4, 5]},
     *     {(3,3), 4, groups: [4, 5, 6]},
     *   ],
     *   groups: [
     *     {[(0,0), (2,0)], numbers: [0], max: 2},
     *     {[(0,1), (1,1), (2,1)], numbers: [0, 1], max: 2},
     *     {[(0,2), (0,3)], numbers: [1, 2], max: 2},
     *     {[(0,4), (1,4)], numbers: [2], max: 2},
     *     {[(2,2), (2,3)], numbers: [1, 2, 3], max: 2},
     *     {[(2,4)], numbers: [2, 3], max: 2},
     *     {[(3,2), (3,4), (4,2), (4,3), (4,4)], numbers: [3], max: 4},
     *   ],
     * }
     *
     * There are 7 groups, meaning combinations will have 7 slots to fill.
     * Combinations: [
     *   (0, 2, 0, 0, 1, 1, 2),
     *   (0, 2, 0, 1, 1, 0, 3),
     *   (0, 2, 1, 0, 0, 1, 3),   #  <- this one shown above
     *   (0, 2, 1, 1, 0, 0, 4),
     *   (1, 1, 0, 0, 2, 0, 2),
     *   (1, 1, 1, 0, 1, 0, 3),
     *   (1, 1, 2, 0, 0, 0, 4),
     * ]
     *
     * These are actually just non-negative integer solutions to the following
     * simultaneous equations:
     *  1) g0 + g1 = 2
     *  2) g1 + g2 + g4 = 3
     *  3) g2 + g3 + g4 + g5 = 2
     *  4) g4 + g5 + g6 = 4
     * subject to:
     *  0 <= g0 <= 2
     *  0 <= g1 <= 2
     *  0 <= g2 <= 2
     *  0 <= g3 <= 2
     *  0 <= g4 <= 2
     *  0 <= g5 <= 2
     *  0 <= g6 <= 4
     *
     * This can alternatively be written as the matrix equation:
     *  |1, 1, 0, 0, 0, 0, 0|        |2|
     *  |0, 1, 1, 0, 1, 0, 0|  . g = |3|
     *  |0, 0, 1, 1, 1, 1, 0|    ~   |2|
     *  |0, 0, 0, 0, 1, 1, 1|        |4|
     *
     * Reducing these equations down gives:
     * (2) - (3): g1 = g3 + g5 + 1  =>  g1 >= 1  and  g3, g5 <= 1
     * (4) - (3): g6 = g2 + g3 + 2  =>  g6 >= 2
     */
}
```
