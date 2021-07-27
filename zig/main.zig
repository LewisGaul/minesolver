//
// Take an input minesweeper board and output info on the solving calculations.
//
// Input format:
//  - Each non-blank line corresponds to a row
//  - Cells in a row are separated by any number of spaces or tabs
//  - A hash symbol (#) represents an unclicked cell
//  - Numbers represent numbers shown on the board
//  - A dot (.) can be used in place of a number 0
//  - An asterisk (*) represents a mine, and may optionally be followed by a
//    number to indicate the number of mines (1 assumed otherwise)
//

const std = @import("std");

const clap = @import("clap");

const ArrayList = std.ArrayList;
const File = std.fs.File;
const Allocator = std.mem.Allocator;
const assert = std.debug.assert;

// -----------------------------------------------------------------------------
// Globals
// -----------------------------------------------------------------------------

pub const log_level = .info;

var allocator: *Allocator = std.heap.page_allocator;

const stdout = std.io.getStdOut().writer();
const stderr = std.io.getStdErr().writer();

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

/// Greatest common divisor.
fn gcd(val1: usize, val2: usize) usize {
    const larger = if (val1 >= val2) val1 else val2;
    const smaller = if (val1 >= val2) val2 else val1;
    return if (larger % smaller == 0)
        smaller
    else
        gcd(smaller, larger % smaller);
}

// Lowest common multiple.
fn lcm(val1: usize, val2: usize) usize {
    if (val1 == val2) return val1;
    if (val1 % val2 == 0) return val1;
    if (val2 % val1 == 0) return val2;
    // TODO: Dodgy int casts...
    return @intCast(u16, val1) * @intCast(u16, val2) / gcd(val1, val2);
}

fn absDifference(val1: anytype, val2: @TypeOf(val1)) @TypeOf(val1) {
    return if (val1 >= val2) val1 - val2 else val2 - val1;
}

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------

const Args = struct {
    input_file: File,
    mines: u16,
    per_cell: u8 = 1,
};

const CellContents = union(enum) {
    // Unclicked cell.
    Unclicked,
    // Revealed space (effectively number 0).
    Space,
    // Revealed number (1, 2, 3, ...).
    Number: u8,
    // Mine or flag (1, 2, 3, ...).
    Mine: u8,

    pub fn format(
        self: @This(),
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        _ = fmt;
        _ = options;
        switch (self) {
            .Unclicked => try writer.writeByte('#'),
            .Space => try writer.writeByte('.'),
            .Number => |n| try writer.print("{}", .{n}),
            .Mine => |n| {
                if (n == 1) {
                    try writer.writeByte('*');
                } else {
                    try writer.print("*{}", .{n});
                }
            },
        }
    }
};

/// Generic rectangular grid type.
///
/// Indexed as follows:
///   0 1 2 3
///   4 5 6 7
fn Grid(comptime T: type) type {
    return struct {
        /// Non-empty square slice of slices.
        data: [][]T,

        const Self = @This();
        const Entry = struct { x: usize, y: usize, value: T };

        /// Iterator for iterating over cells in index order.
        const Iterator = struct {
            grid: Self,
            idx: usize = 0,

            pub fn next(it: *@This()) ?Entry {
                if (it.idx >= it.grid.xSize() * it.grid.ySize())
                    return null;
                const x = it.idx % it.grid.xSize();
                const y = @divFloor(it.idx, it.grid.xSize());
                const entry = Entry{
                    .x = x,
                    .y = y,
                    .value = it.grid.getCell(x, y),
                };
                it.idx += 1;
                return entry;
            }
        };

        /// Allocates memory to copy the data, free by calling 'Grid.deinit()'.
        pub fn fromFlatSlice(x_size: usize, y_size: usize, cells: []const T) !Self {
            if (cells.len != x_size * y_size or cells.len == 0)
                return error.InvalidNumberOfCells;
            const rows: [][]T = try allocator.alloc([]T, y_size);
            for (rows) |*row, y| {
                row.* = try allocator.dupe(T, cells[y * x_size .. (y + 1) * x_size]);
            }
            return Self{ .data = rows };
        }

        pub fn fromOwnedData(data: []const []const T) !Self {
            const new_rows = try allocator.alloc([]T, data.len);
            for (data) |row, i| {
                new_rows[i] = try allocator.dupe(T, row);
            }
            return Self{ .data = new_rows };
        }

        pub fn deinit(g: Self) void {
            for (g.data) |row| {
                allocator.free(row);
            }
            allocator.free(g.data);
        }

        pub fn xSize(g: Self) usize {
            return g.data[0].len;
        }

        pub fn ySize(g: Self) usize {
            return g.data.len;
        }

        /// The 'x' and 'y' arguments must be in range of the grid bounds, as
        /// given by 'Grid.xSize()' and 'Grid.ySize()'.
        pub fn getCell(g: Self, x: usize, y: usize) T {
            return g.data[y][x];
        }

        /// The 'idx' argument must be in range of the grid bounds, as given by
        /// 'Grid.ySize()'.
        pub fn getRow(g: Self, idx: usize) []const T {
            return g.data[idx];
        }

        /// The 'idx' argument must be in range of the grid bounds, as given by
        /// 'Grid.xSize()'.
        pub fn getColumn(g: Self, idx: usize) ![]const T {
            var list = ArrayList(T).init(allocator);
            defer list.deinit();
            try list.ensureTotalCapacity(g.xSize());
            for (g.data) |row| {
                list.appendAssumeCapacity(row[idx]);
            }
            return list.toOwnedSlice();
        }

        pub fn count(g: Self, item: anytype) usize {
            var c: usize = 0;
            for (g.data) |row| {
                for (row) |cell| {
                    if (cell == item) c += 1;
                }
            }
            return c;
        }

        pub fn toStr(g: Self, opts: struct { sep_idx: ?usize = null }) ![]const u8 {
            // Find the max cell width.
            var max_width: u64 = 0;
            var iter = g.iterator();
            while (iter.next()) |entry| {
                const cell_width = std.fmt.count("{}", .{entry.value});
                if (cell_width > max_width) max_width = cell_width;
            }

            var buf = ArrayList(u8).init(allocator);
            defer buf.deinit();
            var writer = buf.writer();
            for (g.data) |row, y| {
                if (y > 0) try writer.writeByte('\n');
                for (row) |cell, i| {
                    if (i > 0) try writer.writeByte(' ');
                    if (opts.sep_idx != null and opts.sep_idx.? == i)
                        try writer.print("| ", .{});
                    const cell_width = std.fmt.count("{}", .{cell});
                    try writer.writeByteNTimes(' ', max_width - cell_width);
                    try writer.print("{}", .{cell});
                }
            }
            return buf.toOwnedSlice();
        }

        /// Make a shallow copy of the grid - memory owned by the caller.
        pub fn copy(g: Self) !Self {
            return Self.fromOwnedData(g.data);
        }

        pub fn iterator(g: Self) Iterator {
            return .{ .grid = g };
        }

        /// Allocates the slice - memory owned by the caller.
        pub fn getNeighbours(g: Self, x: usize, y: usize) ![]const Entry {
            var x_min: usize = x;
            var x_max: usize = x;
            var y_min: usize = y;
            var y_max: usize = y;
            if (x > 0) x_min -= 1;
            if (x < g.xSize() - 1) x_max += 1;
            if (y > 0) y_min -= 1;
            if (y < g.ySize() - 1) y_max += 1;
            const num_nbrs = (x_max - x_min + 1) * (y_max - y_min + 1) - 1;
            var nbrs = ArrayList(Entry).init(allocator);
            try nbrs.ensureTotalCapacity(num_nbrs);
            var i: usize = x_min;
            var j: usize = y_min;
            while (i <= x_max) : (i += 1) {
                while (j <= y_max) : (j += 1) {
                    if (i == x and j == y) continue;
                    nbrs.appendAssumeCapacity(
                        Entry{ .x = i, .y = j, .value = g.getCell(i, j) },
                    );
                }
            }
            return nbrs.toOwnedSlice();
        }
    };
}

const Matrix = struct {
    grid: Grid(isize),

    const Self = @This();

    pub fn init(cells: [][]isize) !Self {
        return Self{ .grid = try Grid(isize).fromOwnedData(cells) };
    }

    pub fn fromFlatSlice(
        comptime T: type,
        x_size: usize,
        y_size: usize,
        cells: []const T,
    ) !Self {
        if (cells.len != x_size * y_size or cells.len == 0)
            return error.InvalidNumberOfCells;
        // TODO: Failing to free memory on error.
        const rows = try allocator.alloc([]isize, y_size);
        for (rows) |*row, y| {
            row.* = try allocator.alloc(isize, x_size);
            for (row.*) |*cell, x| {
                cell.* = cells[y * x_size + x];
            }
        }
        return Self{
            .grid = Grid(isize){ .data = rows },
        };
    }

    /// Return a new matrix containing selected columns from the current matrix.
    pub fn selectColumns(self: Self, col_idxs: []const usize) !Self {
        // TODO: Free memory on error.
        const rows = try allocator.alloc([]isize, self.ySize());
        for (rows) |*row, y| {
            row.* = try allocator.alloc(isize, col_idxs.len);
            for (row.*) |*cell, x| {
                cell.* = self.getCell(col_idxs[x], y);
            }
        }
        return Self{ .grid = Grid(isize){ .data = rows } };
    }

    pub fn fromStr(input: []const u8) !Self {
        var rows = ArrayList([]isize).init(allocator);
        defer rows.deinit();

        var first_line_cols: ?usize = null;
        var lines_iter = std.mem.tokenize(input, "\r\n");
        while (lines_iter.next()) |line| {
            var row = ArrayList(isize).init(allocator);
            defer row.deinit();
            var cells_iter = std.mem.tokenize(line, " \t");
            while (cells_iter.next()) |cell| {
                const val = try std.fmt.parseInt(isize, cell, 10);
                try row.append(val);
            }
            if (row.items.len == 0) {
                // Ignore blank lines.
                continue;
            } else if (first_line_cols == null) {
                first_line_cols = row.items.len;
            } else if (row.items.len != first_line_cols) {
                return error.InconsistentGridShape;
            }
            try rows.append(row.toOwnedSlice());
        }
        if (first_line_cols == null) return error.EmptyInput;

        return Self{ .grid = Grid(isize){ .data = rows.toOwnedSlice() } };
    }

    pub fn deinit(self: Self) void {
        self.grid.deinit();
    }

    pub fn xSize(self: Self) usize {
        return self.grid.xSize();
    }

    pub fn ySize(self: Self) usize {
        return self.grid.ySize();
    }

    /// The 'x' and 'y' arguments must be in range of the grid bounds, as
    /// given by 'xSize()' and 'ySize()'.
    pub fn getCell(self: Self, x: usize, y: usize) isize {
        return self.grid.getCell(x, y);
    }

    /// The 'idx' argument must be in range of the grid bounds, as given by
    /// 'ySize()'.
    pub fn getRow(self: Self, idx: usize) []const isize {
        return self.grid.getRow(idx);
    }

    /// The 'idx' argument must be in range of the grid bounds, as given by
    /// 'xSize()'.
    pub fn getColumn(self: Self, idx: usize) ![]const isize {
        return self.grid.getColumn(idx);
    }

    pub fn toStr(self: Self, opts: struct { sep_idx: ?usize = null }) ![]const u8 {
        return self.grid.toStr(.{ .sep_idx = opts.sep_idx });
    }

    /// Returns a copy of the matrix - memory owned by the caller.
    pub fn copy(self: Self) !Self {
        return Self{ .grid = try self.grid.copy() };
    }

    pub fn matrixAdd(self: Self, other: Self) !Self {
        assert(self.xSize() == other.xSize() and self.ySize() == other.ySize());
        var new_matrix = try self.copy();
        for (new_matrix.grid.data) |row, y| {
            for (row) |*cell, x| {
                cell.* += other.getCell(x, y);
            }
        }
        return new_matrix;
    }

    pub fn matrixSubtract(self: Self, other: Self) !Self {
        assert(self.xSize() == other.xSize() and self.ySize() == other.ySize());
        var new_matrix = try self.copy();
        for (new_matrix.grid.data) |row, y| {
            for (row) |*cell, x| {
                cell.* -= other.getCell(x, y);
            }
        }
        return new_matrix;
    }

    pub fn matrixMultiply(self: Self, other: Self) !Self {
        assert(self.xSize() == other.ySize());
        const new_rows = try allocator.alloc([]isize, self.ySize());
        var y: u8 = 0;
        for (new_rows) |*new_row| {
            new_row.* = try allocator.alloc(isize, other.xSize());
            var x: u8 = 0;
            for (new_row.*) |*cell| {
                cell.* = 0;
                var idx: u8 = 0;
                while (idx < self.xSize()) : (idx += 1) {
                    cell.* += self.getCell(idx, y) * other.getCell(x, idx);
                }
                x += 1;
            }
            y += 1;
        }
        return Self{ .grid = Grid(isize){ .data = new_rows } };
    }

    /// Convert in-place to Reduced-Row-Echelon-Form.
    pub fn rref(self: *Self) void {
        var row_idx: usize = 0;
        var col_idx: usize = 0;
        while (col_idx < self.xSize()) : (col_idx += 1) {
            // std.log.debug("Column {d}, row {d}:\n{s}", .{ col_idx, row_idx, self.toStr(.{}) catch unreachable });
            if (self.findNonZeroRowInColumn(col_idx, row_idx)) |non_zero_row_idx| {
                if (row_idx != non_zero_row_idx) {
                    std.log.debug("Swapping rows {d} and {d}", .{ row_idx, non_zero_row_idx });
                    self.swapRows(row_idx, non_zero_row_idx);
                }
            } else {
                std.log.debug("No row with non-zero value found, skipping column {d}", .{col_idx});
                continue;
            }

            // Now the pivot row 'row_idx' should contain a non-zero value in
            // the pivot column 'col_idx'.
            if (self.getCell(col_idx, row_idx) < 0) {
                std.log.debug("Negating row {d}", .{row_idx});
                self.negateRow(row_idx);
            }
            assert(self.getCell(col_idx, row_idx) > 0);
            const pivot_val = @intCast(usize, self.getCell(col_idx, row_idx));

            var inner_row_idx: usize = 0;
            while (inner_row_idx < self.ySize()) : (inner_row_idx += 1) {
                if (row_idx == inner_row_idx) continue;
                if (self.getCell(col_idx, inner_row_idx) == 0) continue;

                if (self.getCell(col_idx, inner_row_idx) < 0 and inner_row_idx > row_idx) {
                    std.log.debug("Negating row {d}", .{inner_row_idx});
                    self.negateRow(inner_row_idx);
                    assert(self.getCell(col_idx, inner_row_idx) > 0);
                }
                const abs_inner_row_val = std.math.absCast(self.getCell(col_idx, inner_row_idx));

                const pivot_lcm = lcm(pivot_val, abs_inner_row_val);
                const pivot_multiple = pivot_lcm / pivot_val;
                const inner_multiple = pivot_lcm / abs_inner_row_val;
                if (pivot_multiple != 1) {
                    std.log.debug("Multiplying row {d} by {d}", .{ row_idx, pivot_multiple });
                    self.multiplyRow(row_idx, pivot_multiple);
                }
                if (inner_multiple != 1) {
                    std.log.debug("Multiplying row {d} by {d}", .{ inner_row_idx, inner_multiple });
                    self.multiplyRow(inner_row_idx, inner_multiple);
                }
                if (self.getCell(col_idx, inner_row_idx) < 0) {
                    std.log.debug("Adding row {d} to row {d}", .{ row_idx, inner_row_idx });
                    self.addRow(inner_row_idx, row_idx);
                } else {
                    std.log.debug("Subtracting row {d} from row {d}", .{ row_idx, inner_row_idx });
                    self.subtractRow(inner_row_idx, row_idx);
                }

                assert(self.getCell(col_idx, inner_row_idx) == 0);
            }
            row_idx += 1;
        }
    }

    /// Remove all rows following the first row containing only zeros, i.e. all
    /// zero-rows for a matrix in RREF.
    pub fn removeTrailingZeroRows(self: *Self) void {
        var zero_row_idx: ?usize = null;
        for (self.grid.data) |row, y| {
            if (std.mem.allEqual(isize, row, 0)) {
                zero_row_idx = y;
                break;
            }
        }
        if (zero_row_idx) |idx| {
            var y: usize = idx;
            while (y < self.ySize()) : (y += 1) {
                allocator.free(self.grid.data[y]);
            }
            self.grid.data = allocator.shrink(self.grid.data, idx);
        }
    }

    /// Find the first non-zero value in the specified column, returning the row
    /// index.
    fn findNonZeroRowInColumn(self: Self, col_idx: usize, start_row_idx: usize) ?usize {
        var row_idx: usize = start_row_idx;
        while (row_idx < self.ySize()) : (row_idx += 1) {
            if (self.getCell(col_idx, row_idx) != 0) return row_idx;
        }
        return null;
    }

    /// Swap the two rows with the given indexes.
    fn swapRows(self: *Self, row1_idx: usize, row2_idx: usize) void {
        const row1 = self.grid.data[row1_idx];
        const row2 = self.grid.data[row2_idx];
        self.grid.data[row1_idx] = row2;
        self.grid.data[row2_idx] = row1;
    }

    fn negateRow(self: *Self, row_idx: usize) void {
        for (self.grid.data[row_idx]) |*cell| {
            cell.* *= -1;
        }
    }

    fn multiplyRow(self: *Self, row_idx: usize, factor: usize) void {
        if (factor == 1) return;
        for (self.grid.data[row_idx]) |*cell| {
            // TODO: Dodgy int casts...
            cell.* = @intCast(u16, cell.*) * @intCast(u16, factor);
        }
    }

    fn addRow(self: *Self, row_idx: usize, add_row_idx: usize) void {
        for (self.grid.data[row_idx]) |*cell, col_idx| {
            cell.* += self.getCell(col_idx, add_row_idx);
        }
    }

    fn subtractRow(self: *Self, row_idx: usize, sub_row_idx: usize) void {
        for (self.grid.data[row_idx]) |*cell, col_idx| {
            cell.* -= self.getCell(col_idx, sub_row_idx);
        }
    }
};

const Board = struct {
    grid: Grid(CellContents),
    mines: u16,

    const Self = @This();

    /// Allocates memory, but also reuses the provided memory storing the
    /// 'cells' slice, which must all be managed by caller.
    pub fn fromFlatSlice(x_size: usize, y_size: usize, cells: []CellContents, mines: u16) !Self {
        return Self{
            .grid = try Grid(CellContents).fromFlatSlice(x_size, y_size, cells),
            .mines = mines,
        };
    }

    pub fn deinit(self: Self) void {
        self.grid.deinit();
    }

    pub fn xSize(self: Self) usize {
        return self.grid.xSize();
    }

    pub fn ySize(self: Self) usize {
        return self.grid.ySize();
    }

    pub fn toStr(self: Self) ![]const u8 {
        return self.grid.toStr(.{});
    }

    /// Convert a board into a set of simultaneous equations, represented in matrix
    /// form. The memory is owned by the caller.
    pub fn toMatrix(self: Self) !Matrix {
        // Each row of the matrix corresponds to one of the simultaneous equations,
        // which come from each visible number on the board.
        // There is an additional final row corresponding to the equation for the
        // total number of mines in the board.
        // Each column corresponds to an unclicked cell, where the value in the
        // matrix is '1' if the unclicked cell is a neighbour of the row's displayed
        // number, and '0' otherwise.
        // There is a single column on the right corresponding to the RHS of the
        // simultaneous equations, i.e. the value of the number shown in the cell
        // corresponding to that row.
        const num_columns = self.grid.count(.Unclicked) + 1;
        const num_rows = self.grid.count(.Number) + 1;
        const all_matrix_cells: []isize = try allocator.alloc(isize, num_columns * num_rows);
        defer allocator.free(all_matrix_cells);

        var j: usize = 0; // Matrix row index
        var iter1 = self.grid.iterator();
        // Iterate over the numbers in the board, which correspond to the rows of
        // the matrix.
        while (iter1.next()) |num_entry| {
            if (num_entry.value != .Number) continue;
            const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
            var i: usize = 0; // Matrix column index
            var iter2 = self.grid.iterator();
            var nbr_mines: usize = 0;
            while (iter2.next()) |nbr_entry| {
                const is_nbr_x = absDifference(nbr_entry.x, num_entry.x) <= 1;
                const is_nbr_y = absDifference(nbr_entry.y, num_entry.y) <= 1;
                const is_nbr = is_nbr_x and is_nbr_y;
                switch (nbr_entry.value) {
                    .Mine => |m| {
                        if (is_nbr) nbr_mines += m;
                        continue;
                    },
                    .Unclicked => {}, // Fall through
                    else => continue,
                }
                row[i] = if (is_nbr) 1 else 0;
                i += 1;
            }
            assert(i == num_columns - 1);
            // Add the final column value - the RHS of the equation.
            if (num_entry.value.Number < nbr_mines) return error.TooManyMinesAroundNumber;
            row[i] = num_entry.value.Number - @intCast(isize, nbr_mines);
            j += 1;
        }
        assert(j == num_rows - 1);

        // Add the final row on the end, corresponding to the total number of mines.
        var i: usize = 0;
        const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
        while (i < num_columns - 1) : (i += 1) {
            row[i] = 1;
        }
        row[i] = self.mines;

        return Matrix.fromFlatSlice(isize, num_columns, num_rows, all_matrix_cells);
    }
};

const RectangularIterator = struct {
    max_vals: []u16,
    idx: usize = 0,

    pub fn next(self: *@This()) !?[]const u16 {
        if (self.idx >= self.size()) return null;
        defer self.idx += 1;
        const result = try allocator.alloc(u16, self.max_vals.len);
        var counter: usize = 1;
        for (self.max_vals) |val, i| {
            result[i] = @intCast(u16, (self.idx / counter) % (val + 1));
            counter *= val + 1;
        }
        return result;
    }

    pub fn size(self: @This()) usize {
        var result: usize = 1;
        for (self.max_vals) |val| {
            result *= val + 1;
        }
        return result;
    }
};

const Solver = struct {
    board: Board,
    per_cell: u8,
    /// Each inner slice contains indexes to grid positions.
    groups: []const []const usize,
    matrix: Matrix,

    const Self = @This();

    /// The board is still owned by the caller and should be independently
    /// deinitialised.
    pub fn init(board: Board, per_cell: u8) !Self {
        var self = Self{
            .board = board,
            .per_cell = per_cell,
            .groups = undefined,
            .matrix = try board.toMatrix(),
        };
        try self.reduceToGroups();
        self.matrix.rref();
        self.matrix.removeTrailingZeroRows();
        return self;
    }

    pub fn deinit(self: Self) void {
        for (self.groups) |group| {
            allocator.free(group);
        }
        allocator.free(self.groups);
        self.matrix.deinit();
    }

    /// Returns a slice containing valid mine configurations, which are
    /// represented as slices of length equal to the number of groups (i.e. the
    /// number of matrix columns), with the values being the number of mines in
    /// the corresponding group.
    ///
    /// Memory owned by the caller.
    pub fn solve(self: Self) ![]const []const u16 {
        // We begin with a matrix in RREF, such as:
        //   1  0  0  0  0  0  1  1 |  1
        //   0  1  0  0  0  1  0  0 |  4
        //   0  0  1  0  0  0 -1 -1 |  1
        //   0  0  0  1  0 -1  1  0 | -2
        //   0  0  0  0  1  1  0  1 |  4
        //
        // We need to transform this to something like:
        //   | g1 |   |  1 |   |  0  1  1 |
        //   | g2 |   |  4 |   |  1  0  0 |   | g6 |
        //   | g3 | = |  1 | - |  0 -1 -1 | . | g7 |
        //   | g4 |   | -2 |   | -1  1  0 |   | g8 |
        //   | g5 |   |  4 |   |  1  0  1 |
        //
        // We then iterate over possible values for the free variables (g6, g7,
        // g8) and matrix multiply to find the fixed variables, checking they're
        // in range.
        //

        // Check for inconsistent matrix - row with all-zero LHS, non-zero RHS
        // (last row in RREF with zero-rows removed).
        {
            const last_row = self.matrix.grid.data[self.matrix.ySize() - 1];
            if (std.mem.allEqual(isize, last_row[0 .. self.matrix.xSize() - 1], 0) and
                last_row[self.matrix.xSize() - 1] != 0)
            {
                return error.InvalidMatrixEquations;
            }
        }

        const rhs_vec = try self.matrix.selectColumns(&[_]usize{self.matrix.xSize() - 1});
        const col_categorisation = try self.categoriseColumns();
        const fixed_col_idxs = col_categorisation.fixed;
        const free_col_idxs = col_categorisation.free;

        var configs = ArrayList([]u16).init(allocator);
        defer configs.deinit();

        if (free_col_idxs.len == 0) {
            const config = try allocator.alloc(u16, self.groups.len);
            for (config) |*val, i| {
                val.* = @intCast(u16, rhs_vec.getCell(0, i));
            }
            try configs.append(config);
            return configs.toOwnedSlice();
        }

        const free_var_matrix = try self.matrix.selectColumns(free_col_idxs);

        const max_free_vals = try allocator.alloc(u16, free_col_idxs.len);
        // TODO: Improve this - big performance sink!
        // To start with, just iterate over the full range for each group, i.e.
        // 0 to the size of the group multiplied by max per cell.
        for (max_free_vals) |*val, i| {
            val.* = @intCast(u16, self.groups[free_col_idxs[i]].len) * self.per_cell;
        }
        var iter = RectangularIterator{ .max_vals = max_free_vals };
        while (try iter.next()) |free_vals| {
            defer allocator.free(free_vals);
            const free_var_vec = try Matrix.fromFlatSlice(
                u16,
                1,
                free_vals.len,
                free_vals,
            );
            defer free_var_vec.deinit();

            const fixed_var_vec = blk: {
                const multiplied_matrix = try free_var_matrix.matrixMultiply(free_var_vec);
                defer multiplied_matrix.deinit();
                break :blk try rhs_vec.matrixSubtract(multiplied_matrix);
            };
            defer fixed_var_vec.deinit();

            var invalid = false;
            for (fixed_var_vec.grid.data) |row, idx| {
                if (row[0] < 0 or
                    row[0] > self.groups[fixed_col_idxs[idx]].len * self.per_cell)
                {
                    invalid = true;
                    break;
                }
            }
            if (invalid) continue;

            var config = try allocator.alloc(u16, self.groups.len);
            errdefer allocator.free(config);

            for (fixed_col_idxs) |grp_idx, fixed_col_idx| {
                config[grp_idx] = @intCast(u16, fixed_var_vec.getCell(0, fixed_col_idx));
            }
            for (free_col_idxs) |grp_idx, free_col_idx| {
                config[grp_idx] = @intCast(u16, free_var_vec.getCell(0, free_col_idx));
            }
            try configs.append(config);
        }

        if (configs.items.len == 0) return error.InvalidMatrixEquations;

        return configs.toOwnedSlice();
    }

    fn reduceToGroups(self: *Self) !void {
        var groups = ArrayList([]const usize).init(allocator);
        defer groups.deinit();

        // Indexes of columns to remove.
        var remove_columns = ArrayList(usize).init(allocator);
        defer remove_columns.deinit();

        // Iterate over columns to find groups.
        var x1: usize = 0;
        while (x1 < self.matrix.xSize() - 2) : (x1 += 1) {
            // If already included in a group, skip over this column.
            if (std.mem.indexOfScalar(usize, remove_columns.items, x1)) |_| continue;

            const column1 = try self.matrix.getColumn(x1);
            defer allocator.free(column1);

            var group = ArrayList(usize).init(allocator);
            defer group.deinit();
            try group.append(x1);

            var x2: usize = x1 + 1;
            while (x2 < self.matrix.xSize() - 1) : (x2 += 1) {
                const column2 = try self.matrix.getColumn(x2);
                defer allocator.free(column2);
                if (std.mem.eql(isize, column1, column2)) {
                    try remove_columns.append(x2);
                    try group.append(x2);
                }
            }
            try groups.append(group.toOwnedSlice());
        }

        const new_x_size = self.matrix.xSize() - remove_columns.items.len;
        const y_size = self.matrix.ySize();

        var new_matrix_cells = ArrayList(isize).init(allocator);
        defer new_matrix_cells.deinit();
        try new_matrix_cells.ensureTotalCapacity(new_x_size * y_size);

        var iter = self.matrix.grid.iterator();
        while (iter.next()) |entry| {
            if (std.mem.indexOfScalar(usize, remove_columns.items, entry.x)) |_| continue;
            new_matrix_cells.appendAssumeCapacity(entry.value);
        }

        self.groups = groups.toOwnedSlice();
        self.matrix.deinit();
        self.matrix = try Matrix.fromFlatSlice(
            isize,
            new_x_size,
            y_size,
            new_matrix_cells.items,
        );
    }

    const ColumnCategorisation = struct { fixed: []const usize, free: []const usize };

    fn categoriseColumns(self: Self) !ColumnCategorisation {
        var fixed_cols = ArrayList(usize).init(allocator);
        defer fixed_cols.deinit();
        var free_cols = ArrayList(usize).init(allocator);
        defer free_cols.deinit();

        var x: usize = 0;
        var y: usize = 0;
        while (y < self.matrix.ySize()) : (y += 1) {
            while (x < self.matrix.xSize() - 1) : (x += 1) {
                if (self.matrix.getCell(x, y) == 0) {
                    try free_cols.append(x);
                } else {
                    try fixed_cols.append(x);
                    x += 1;
                    break;
                }
            }
        }
        // All remaining columns are free columns (leading '1's all found).
        while (x < self.matrix.xSize() - 1) : (x += 1) {
            try free_cols.append(x);
        }
        return ColumnCategorisation{
            .fixed = fixed_cols.toOwnedSlice(),
            .free = free_cols.toOwnedSlice(),
        };
    }
};

// -----------------------------------------------------------------------------
// Input parsing
// -----------------------------------------------------------------------------

fn parseInputCell(input: []const u8) !CellContents {
    switch (input[0]) {
        '1'...'9' => return CellContents{
            .Number = try std.fmt.parseUnsigned(u8, input, 10),
        },
        '0', '.' => {
            if (input.len > 1) return error.UnexpectedCellText;
            return CellContents.Space;
        },
        '#' => {
            if (input.len > 1) return error.UnexpectedCellText;
            return CellContents.Unclicked;
        },
        '*' => {
            if (input.len == 1) return CellContents{ .Mine = 1 };
            if (input[1] == '0') return error.UnexpectedCellText;
            return CellContents{ .Mine = try std.fmt.parseUnsigned(u8, input[1..], 10) };
        },
        else => return error.UnexpectedCellText,
    }
}

fn parseInputBoard(input: []const u8, mines: u16) !Board {
    var cells_array = ArrayList(CellContents).init(allocator);
    defer cells_array.deinit();

    var rows: u8 = 0;
    var first_line_cols: ?u8 = null;
    var lines_iter = std.mem.tokenize(input, "\r\n");
    while (lines_iter.next()) |line| {
        var cells_iter = std.mem.tokenize(line, " \t");
        var cols: u8 = 0;
        while (cells_iter.next()) |cell| {
            try cells_array.append(try parseInputCell(cell));
            cols += 1;
        }
        if (cols == 0) {
            // Ignore blank lines.
            continue;
        } else if (first_line_cols == null) {
            first_line_cols = cols;
        } else if (cols != first_line_cols) {
            return error.InconsistentGridShape;
        }
        rows += 1;
    }
    if (first_line_cols == null) return error.EmptyInput;

    return Board.fromFlatSlice(
        first_line_cols.?,
        rows,
        cells_array.items,
        mines,
    );
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

fn parseArgs() !Args {
    // First we specify what parameters our program can take.
    // We can use 'parseParam()' to parse a string to a 'Param(Help)'.
    const params = comptime [_]clap.Param(clap.Help){
        clap.parseParam("-h, --help           Display this help and exit") catch unreachable,
        clap.parseParam("<MINES>              Number of mines") catch unreachable,
        clap.parseParam("-f, --file <PATH>    Input file (defaults to stdin)") catch unreachable,
        clap.parseParam("--per-cell <NUM>     Max number of mines per cell") catch unreachable,
    };

    // Initalize diagnostics for reporting parsing errors.
    var diag = clap.Diagnostic{};
    var clap_args = clap.parse(clap.Help, &params, .{ .diagnostic = &diag }) catch |err| {
        diag.report(stderr, err) catch {};
        return err;
    };
    defer clap_args.deinit();

    if (clap_args.flag("--help")) {
        try stderr.print("{s} ", .{clap_args.exe_arg});
        try clap.usage(stderr, &params);
        try stderr.writeByte('\n');
        try clap.help(stderr, &params);
        std.process.exit(0);
    }

    if (clap_args.positionals().len != 1) {
        try stderr.writeAll("Expected exactly one positional arg\n");
        return error.InvalidArgument;
    }

    const input_file = blk: {
        if (clap_args.option("--file")) |file| {
            break :blk std.fs.cwd().openFile(file, .{}) catch |err| {
                try stderr.print("Failed to open file {s}\n", .{file});
                return err;
            };
        } else break :blk std.io.getStdIn();
    };

    const mines_arg = clap_args.positionals()[0];
    const mines = std.fmt.parseUnsigned(u16, mines_arg, 10) catch |err| {
        try stderr.writeAll("Expected positive integer number of mines\n");
        return err;
    };
    if (mines == 0) {
        try stderr.writeAll("Expected positive integer number of mines\n");
        return error.InvalidArgument;
    }

    var args = Args{ .input_file = input_file, .mines = mines };

    if (clap_args.option("--per-cell")) |per_cell| {
        args.per_cell = try std.fmt.parseUnsigned(u8, per_cell, 10);
        if (args.per_cell == 0) {
            try stderr.writeAll("Max number of mines per cell must be greater than 0\n");
            return error.InvalidArgument;
        }
    }

    return args;
}

pub fn main() !u8 {
    // Set up an allocator - no need to free memory as we go.
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    allocator = &gpa.allocator;

    const args = parseArgs() catch |err| switch (err) {
        error.InvalidArgument,
        error.MissingValue,
        error.DoesntTakeValue,
        => return 2,
        else => return 1,
    };

    try stdout.print(
        "Parsed args:\nmines={d}, per_cell={d}\n",
        .{ args.mines, args.per_cell },
    );

    const max_size = 1024 * 1024; // 1MB
    const input = args.input_file.readToEndAlloc(
        allocator,
        max_size,
    ) catch |err| switch (err) {
        error.FileTooBig => {
            try stderr.print("Failed to read input, {s} - 1MB max\n", .{@errorName(err)});
            return 1;
        },
        else => {
            try stderr.print("Failed to read input, {s}\n", .{@errorName(err)});
            return 1;
        },
    };

    try stdout.print("\n", .{});

    const board = try parseInputBoard(input, args.mines);
    defer board.deinit();
    try stdout.print("Board:\n{s}\n", .{try board.toStr()});

    try stdout.print("\n", .{});

    var matrix = try board.toMatrix();
    defer matrix.deinit();
    try stdout.print(
        "Matrix:\n{s}\n",
        .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
    );

    try stdout.print("\n", .{});

    matrix.rref();
    try stdout.print(
        "RREF matrix:\n{s}\n",
        .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
    );

    try stdout.print("\n", .{});

    const solver = try Solver.init(board, args.per_cell);
    defer solver.deinit();
    try stdout.print(
        "Solver matrix:\n{s}\n",
        .{try solver.matrix.toStr(.{ .sep_idx = solver.matrix.xSize() - 1 })},
    );

    try stdout.print("\n", .{});

    try stdout.print("Solver groups:\n", .{});
    for (solver.groups) |group, i| {
        try stdout.print("{d}: {d}\n", .{ i, group });
    }

    try stdout.print("\n", .{});

    const configs = try solver.solve();
    try stdout.print("Mine configurations:\n", .{});
    for (configs) |cfg, i| {
        try stdout.print("{d}: {d}\n", .{ i, cfg });
    }

    return 0;
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

test "Grid init/deinit" {
    allocator = std.testing.allocator;
    const data = [_]u8{ 1, 2, 3, 4, 5, 6 };
    const grid = try Grid(u8).fromFlatSlice(2, 3, &data);
    grid.deinit();
}

test "Matrix multiply" {
    allocator = std.testing.allocator;
    const matrix1 = try Matrix.fromStr(
        \\ 1 2 3
        \\ 4 5 6
    );
    defer matrix1.deinit();

    const matrix2 = try Matrix.fromStr(
        \\ 1 2
        \\ 4 5
        \\ 7 8
    );
    defer matrix2.deinit();

    const matrix3 = try matrix1.matrixMultiply(matrix2);
    defer matrix3.deinit();

    const mat_str = try matrix3.toStr(.{});
    defer allocator.free(mat_str);
    const exp_str =
        \\30 36
        \\66 81
    ;
    try std.testing.expectEqualStrings(exp_str, mat_str);
}

test "Matrix remove trailing zero rows" {
    allocator = std.testing.allocator;
    var matrix = try Matrix.fromStr(
        \\ 1 2 3
        \\ 4 5 6
        \\ 0 0 0
        \\ 0 0 0
    );
    defer matrix.deinit();

    matrix.removeTrailingZeroRows();

    const mat_str = try matrix.toStr(.{});
    defer allocator.free(mat_str);
    const exp_str =
        \\1 2 3
        \\4 5 6
    ;
    try std.testing.expectEqualStrings(exp_str, mat_str);
}

test "Rectangular iterator" {
    allocator = std.testing.allocator;
    var max_vals = [_]u16{ 2, 0, 1 };
    var iter = RectangularIterator{ .max_vals = &max_vals };
    try std.testing.expectEqual(@as(usize, 6), iter.size());

    const expected_vals_array = [_][]const u16{
        &[_]u16{ 0, 0, 0 },
        &[_]u16{ 1, 0, 0 },
        &[_]u16{ 2, 0, 0 },
        &[_]u16{ 0, 0, 1 },
        &[_]u16{ 1, 0, 1 },
        &[_]u16{ 2, 0, 1 },
    };
    for (expected_vals_array) |exp_vals| {
        if (try iter.next()) |vals| {
            try std.testing.expectEqualSlices(u16, exp_vals, vals);
            allocator.free(vals);
        } else return error.ExpectedIteratorItem;
    }
    try std.testing.expectEqual(@as(?[]const u16, null), try iter.next());
}

test "Solver: invalid board" {
    allocator = std.testing.allocator;
    const board = try parseInputBoard(
        \\ # 1 2
        \\ # # #
    , 2);
    defer board.deinit();

    const solver = try Solver.init(board, 1);
    defer solver.deinit();

    try std.testing.expectError(error.InvalidMatrixEquations, solver.solve());
}
