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

var allocator: *Allocator = undefined;

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

/// Generic square grid type.
///
/// Indexed as follows:
///   0 1 2 3
///   4 5 6 7
fn Grid(comptime T: type) type {
    return struct {
        /// Non-empty square slice of slices.
        data: [][]T,

        const Self = @This();
        const Entry = struct { x: u8, y: u8, value: T };

        /// Iterator for iterating over cells in index order.
        const Iterator = struct {
            grid: Self,
            idx: u8 = 0,

            pub fn next(it: *@This()) ?Entry {
                if (it.idx >= it.grid.xSize() * it.grid.ySize()) return null;
                const x = it.idx % it.grid.xSize();
                const y = @divFloor(it.idx, it.grid.xSize());
                const entry = Entry{
                    .x = x,
                    .y = y,
                    .value = it.grid.get(x, y),
                };
                it.idx += 1;
                return entry;
            }
        };

        /// Allocates memory, but also reuses the provided memory storing the
        /// 'cells' slice, which must all be managed by caller.
        pub fn fromFlatSlice(x_size: u8, y_size: u8, cells: []T) !Self {
            if (cells.len != x_size * y_size or cells.len == 0) return error.InvalidNumberOfCells;
            const rows: [][]T = try allocator.alloc([]T, y_size);
            var y: u8 = 0;
            while (y < y_size) : (y += 1) {
                rows[y] = cells[y * x_size .. (y + 1) * x_size];
            }
            return Self{ .data = rows };
        }

        pub fn xSize(g: Self) u8 {
            return std.math.cast(u8, g.data[0].len) catch unreachable;
        }

        pub fn ySize(g: Self) u8 {
            return std.math.cast(u8, g.data.len) catch unreachable;
        }

        /// The 'x' and 'y' arguments must be in range of the grid bounds, as
        /// given by 'Grid.xSize()' and 'Grid.ySize()'.
        pub fn get(g: Self, x: u8, y: u8) T {
            return g.data[y][x];
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

        pub fn toStr(g: Self, opts: struct { sep_idx: ?u8 = null }) ![]const u8 {
            // Find the max cell width.
            var max_width: u64 = 0;
            var iter = g.iterator();
            while (iter.next()) |entry| {
                const cell_width = std.fmt.count("{}", .{entry.value});
                if (cell_width > max_width) max_width = cell_width;
            }

            var buf = ArrayList(u8).init(allocator);
            errdefer buf.deinit();
            var writer = buf.writer();
            for (g.data) |row, y| {
                if (y > 0) try writer.writeByte('\n');
                for (row) |cell, i| {
                    if (opts.sep_idx != null and opts.sep_idx.? == i)
                        try writer.print("| ", .{});
                    const cell_width = std.fmt.count("{}", .{cell});
                    try writer.writeByteNTimes(' ', max_width - cell_width);
                    try writer.print("{} ", .{cell});
                }
            }
            return buf.toOwnedSlice();
        }

        pub fn iterator(g: Self) Iterator {
            return .{ .grid = g };
        }

        /// Allocates the slice - memory owned by the caller.
        pub fn getNeighbours(g: Self, x: u8, y: u8) ![]const Entry {
            var x_min: u8 = x;
            var x_max: u8 = x;
            var y_min: u8 = y;
            var y_max: u8 = y;
            if (x > 0) x_min -= 1;
            if (x < g.xSize() - 1) x_max += 1;
            if (y > 0) y_min -= 1;
            if (y < g.ySize() - 1) y_max += 1;
            const num_nbrs = (x_max - x_min + 1) * (y_max - y_min + 1) - 1;
            var nbrs = ArrayList(Entry).init(allocator);
            try nbrs.ensureTotalCapacity(num_nbrs);
            var i: u8 = x_min;
            var j: u8 = y_min;
            while (i <= x_max) : (i += 1) {
                while (j <= y_max) : (j += 1) {
                    if (i == x and j == y) continue;
                    nbrs.appendAssumeCapacity(
                        Entry{ .x = i, .y = j, .value = g.get(i, j) },
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

    pub fn init(cells: [][]isize) Self {
        return .{ .grid = .{ .data = cells } };
    }

    pub fn xSize(self: Self) u8 {
        return self.grid.xSize();
    }

    pub fn ySize(self: Self) u8 {
        return self.grid.ySize();
    }

    /// The 'x' and 'y' arguments must be in range of the grid bounds, as
    /// given by 'xSize()' and 'ySize()'.
    pub fn get(self: Self, x: u8, y: u8) isize {
        return self.grid.get(x, y);
    }

    pub fn toStr(self: Self, opts: struct { sep_idx: ?u8 = null }) ![]const u8 {
        return self.grid.toStr(.{ .sep_idx = opts.sep_idx });
    }

    /// Convert in-place to Reduced-Row-Echelon-Form.
    pub fn rref(self: *Self) void {
        var row_idx: u8 = 0;
        var col_idx: u8 = 0;
        while (col_idx < self.xSize()) : (col_idx += 1) {
            std.log.debug("Column {d}, row {d}:\n{s}", .{ col_idx, row_idx, self.toStr(.{}) catch unreachable });
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
            if (self.get(col_idx, row_idx) < 0) {
                std.log.debug("Negating row {d}", .{row_idx});
                self.negateRow(row_idx);
            }
            assert(self.get(col_idx, row_idx) > 0);
            const pivot_val = @intCast(usize, self.get(col_idx, row_idx));

            var inner_row_idx: u8 = 0;
            while (inner_row_idx < self.ySize()) : (inner_row_idx += 1) {
                if (row_idx == inner_row_idx) continue;
                if (self.get(col_idx, inner_row_idx) == 0) continue;

                if (self.get(col_idx, inner_row_idx) < 0 and inner_row_idx > row_idx) {
                    std.log.debug("Negating row {d}", .{inner_row_idx});
                    self.negateRow(inner_row_idx);
                    assert(self.get(col_idx, inner_row_idx) > 0);
                }
                const abs_inner_row_val = std.math.absCast(self.get(col_idx, inner_row_idx));

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
                if (self.get(col_idx, inner_row_idx) < 0) {
                    std.log.debug("Adding row {d} to row {d}", .{ row_idx, inner_row_idx });
                    self.addRow(inner_row_idx, row_idx);
                } else {
                    std.log.debug("Subtracting row {d} from row {d}", .{ row_idx, inner_row_idx });
                    self.subtractRow(inner_row_idx, row_idx);
                }

                assert(self.get(col_idx, inner_row_idx) == 0);
            }
            row_idx += 1;
        }
    }

    /// Find the first non-zero value in the specified column, returning the row
    /// index.
    fn findNonZeroRowInColumn(self: Self, col_idx: u8, start_row_idx: u8) ?u8 {
        var row_idx: u8 = start_row_idx;
        while (row_idx < self.ySize()) : (row_idx += 1) {
            if (self.get(col_idx, row_idx) != 0) return row_idx;
        }
        return null;
    }

    /// Swap the two rows with the given indexes.
    fn swapRows(self: *Self, row1_idx: u8, row2_idx: u8) void {
        const row1 = self.grid.data[row1_idx];
        const row2 = self.grid.data[row2_idx];
        self.grid.data[row1_idx] = row2;
        self.grid.data[row2_idx] = row1;
    }

    fn negateRow(self: *Self, row_idx: u8) void {
        for (self.grid.data[row_idx]) |*cell| {
            cell.* *= -1;
        }
    }

    fn multiplyRow(self: *Self, row_idx: u8, factor: usize) void {
        if (factor == 1) return;
        for (self.grid.data[row_idx]) |*cell| {
            // TODO: Dodgy int casts...
            cell.* = @intCast(u16, cell.*) * @intCast(u16, factor);
        }
    }

    fn addRow(self: *Self, row_idx: u8, sub_row_idx: u8) void {
        for (self.grid.data[row_idx]) |*cell, col_idx| {
            cell.* += self.get(@intCast(u8, col_idx), sub_row_idx);
        }
    }

    fn subtractRow(self: *Self, row_idx: u8, sub_row_idx: u8) void {
        for (self.grid.data[row_idx]) |*cell, col_idx| {
            cell.* -= self.get(@intCast(u8, col_idx), sub_row_idx);
        }
    }
};

const Board = struct {
    grid: Grid(CellContents),
    mines: u16,

    const Self = @This();

    /// Allocates memory, but also reuses the provided memory storing the
    /// 'cells' slice, which must all be managed by caller.
    pub fn fromFlatSlice(x_size: u8, y_size: u8, cells: []CellContents, mines: u16) !Self {
        return Self{
            .grid = try Grid(CellContents).fromFlatSlice(x_size, y_size, cells),
            .mines = mines,
        };
    }

    pub fn xSize(self: Self) u8 {
        return self.grid.xSize();
    }

    pub fn ySize(self: Self) u8 {
        return self.grid.ySize();
    }

    pub fn toStr(self: Self) ![]const u8 {
        return self.grid.toStr(.{});
    }

    /// Convert a board into a set of simultaneous equations, represented in matrix
    /// form. The memory is owned by the caller.
    pub fn toMatrix(self: Self) !Matrix {
        const absInt = std.math.absInt;
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
        const rows: [][]isize = try allocator.alloc([]isize, num_rows);
        const all_matrix_cells: []isize = try allocator.alloc(isize, num_columns * num_rows);

        var j: u8 = 0; // Matrix row index
        var iter1 = self.grid.iterator();
        // Iterate over the numbers in the board, which correspond to the rows of
        // the matrix.
        while (iter1.next()) |num_entry| {
            if (num_entry.value != .Number) continue;
            const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
            var i: u8 = 0; // Matrix column index
            var iter2 = self.grid.iterator();
            var nbr_mines: u8 = 0;
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
            row[i] = num_entry.value.Number - nbr_mines;
            rows[j] = row;
            j += 1;
        }
        assert(j == num_rows - 1);

        // Add the final row on the end, corresponding to the total number of mines.
        var i: u8 = 0;
        const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
        while (i < num_columns - 1) : (i += 1) {
            row[i] = 1;
        }
        row[i] = self.mines;
        rows[j] = row;

        return Matrix.init(rows);
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
    errdefer cells_array.deinit();

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
        cells_array.toOwnedSlice(),
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
    // Use an arena allocator - no need to free memory as we go.
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    allocator = &arena.allocator;

    const args = parseArgs() catch |err| switch (err) {
        error.InvalidArgument,
        error.MissingValue,
        error.DoesntTakeValue,
        => return 2,
        else => return 1,
    };

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

    const board = try parseInputBoard(input, args.mines);
    std.debug.print("Board:\n{s}\n", .{try board.toStr()});

    std.debug.print("\n", .{});
    var matrix = try board.toMatrix();
    std.debug.print(
        "Matrix:\n{s}\n",
        .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
    );

    std.debug.print("\n", .{});
    matrix.rref();
    std.debug.print(
        "RREF matrix:\n{s}\n",
        .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
    );

    return 0;
}
