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

var start_milli_time: u64 = 0;

// This is the log level read by std.log - without this it's impossible to get
// info/debug logging in release build modes, and we control it with the
// --verbose CLI flag anyway.
pub const log_level: std.log.Level = .debug;

var internal_log_level: std.log.Level = .info;

/// Override the default logger in std.log, including a timestamp in the output.
pub fn log(
    comptime message_level: std.log.Level,
    comptime scope: @Type(.EnumLiteral),
    comptime format: []const u8,
    args: anytype,
) void {
    if (@enumToInt(message_level) <= @enumToInt(internal_log_level)) {
        const milli_time = @intCast(u64, std.time.milliTimestamp()) - start_milli_time;
        const level_txt = switch (message_level) {
            // zig fmt: off
            .emerg  => "EMERG",
            .alert  => "ALERT",
            .crit   => " CRIT",
            .err    => "ERROR",
            .warn   => " WARN",
            .notice => " NOTE",
            .info   => " INFO",
            .debug  => "DEBUG",
            // zig fmt: on
        };
        const prefix2 = if (scope == .default) ": " else "(" ++ @tagName(scope) ++ "): ";
        const held = std.debug.getStderrMutex().acquire();
        defer held.release();
        nosuspend stderr.print("{d:>2}.{d:0>3} ", .{ milli_time / 1000, milli_time % 1000 }) catch return;
        nosuspend stderr.print("[" ++ level_txt ++ "]" ++ prefix2 ++ format ++ "\n", args) catch return;
    }
}

var allocator: *Allocator = std.heap.page_allocator;

// Have to get these at runtime on Windows.
var stdout: std.fs.File.Writer = undefined;
var stderr: std.fs.File.Writer = undefined;

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

/// Returns the log of the number of combinations for group of size 's',
/// 'm' mines, and 'xmax' mines per cell.
fn logCombs(s: usize, m: usize, xmax: u8) f64 {
    assert(m <= s * xmax);
    if (s == 1) return 0;
    if (xmax >= m) {
        // s**m
        return std.math.ln(@intToFloat(f64, s)) * @intToFloat(f64, m);
    } else if (xmax == 1) {
        // Falling factorial s!/m!
        var result: f64 = 0;
        var i: usize = s - m + 1;
        while (i <= s) : (i += 1) {
            result += std.math.ln(@intToFloat(f64, i));
        }
        return result;
    } else if (xmax == 2) {
        // A sum expression... see docs.
        // Using floats accepting loss of precision, since the number of
        // combinations can easily exceed something like u1000 max. Can't do
        // all calculations with logs as we need to add numbers of combinations!
        var tot: f64 = 0;
        var d: usize = if (s >= m) 0 else m - s; // Number of double mines
        while (d <= m / 2) : (d += 1) {
            var ln_val: f64 = 0;
            { // s! / (s - m + d)!
                var i: usize = s + d - m + 1;
                while (i <= s) : (i += 1) {
                    ln_val += std.math.ln(@intToFloat(f64, i));
                }
            }
            { // m! / (m - 2*d)!
                var i: usize = m - 2 * d + 1;
                while (i <= m) : (i += 1) {
                    ln_val += std.math.ln(@intToFloat(f64, i));
                }
            }
            { // 1 / (2**d * d!)
                ln_val -= std.math.ln(2.0) * @intToFloat(f64, d);
                var i: usize = 1;
                while (i <= d) : (i += 1) {
                    ln_val -= std.math.ln(@intToFloat(f64, i));
                }
            }
            tot += std.math.exp(ln_val);
        }
        return std.math.ln(tot);
    } else if (xmax == 3) {
        // A horrible nested sum expression... see docs.
        // Using floats accepting loss of precision, since the number of
        // combinations can easily exceed something like u1000 max. Can't do
        // all calculations with logs as we need to add numbers of combinations!
        var tot: f64 = 0;
        var t: usize = if (2 * s >= m) 0 else m - 2 * s; // Number of triple mines
        while (t <= m / 3) : (t += 1) {
            var d: usize = if (s >= m - 2 * t) 0 else m - 2 * t - s; // Number of double mines
            while (d <= (m - 3 * t) / 2) : (d += 1) {
                var ln_val: f64 = 0;
                { // s! / (s - m + d + 2*t)!
                    var i: usize = s + d + 2 * t - m + 1;
                    while (i <= s) : (i += 1) {
                        ln_val += std.math.ln(@intToFloat(f64, i));
                    }
                }
                { // m! / (m - 2*d - 3*t)!
                    var i: usize = m - 2 * d - 3 * t + 1;
                    while (i <= m) : (i += 1) {
                        ln_val += std.math.ln(@intToFloat(f64, i));
                    }
                }
                { // 1 / ( (2!)**d * (3!)**t )
                    ln_val -= std.math.ln(2.0) * @intToFloat(f64, d);
                    ln_val -= std.math.ln(6.0) * @intToFloat(f64, t);
                }
                { // 1/d!
                    var i: usize = 1;
                    while (i <= d) : (i += 1) {
                        ln_val -= std.math.ln(@intToFloat(f64, i));
                    }
                }
                { // 1/t!
                    var i: usize = 1;
                    while (i <= t) : (i += 1) {
                        ln_val -= std.math.ln(@intToFloat(f64, i));
                    }
                }
                tot += std.math.exp(ln_val);
            }
        }
        return std.math.ln(tot);
    } else {
        @panic("Unable to calculate number of combinations for per_cell > 3");
    }
    unreachable;
}

/// Calculate the probability a cell contains a mine in a group of size 's'
/// containing 'm' mines and with max per cell of 'xmax'.
fn unsafeProb(s: usize, m: usize, xmax: u8) f64 {
    if (m > s * xmax) return 0;
    if (xmax == 1) {
        return @intToFloat(f64, m) / @intToFloat(f64, s);
    } else if (xmax >= m) {
        return 1 - std.math.pow(f64, 1 - 1.0 / @intToFloat(f64, s), @intToFloat(f64, m));
    } else if (m > xmax * (s - 1)) {
        return 1;
    } else {
        return 1 - std.math.exp(logCombs(s - 1, m, xmax) - logCombs(s, m, xmax));
    }
    unreachable;
}

// -----------------------------------------------------------------------------
// Structs
// -----------------------------------------------------------------------------

const MinesInfo = union(enum) {
    Num: u16,
    Density: f64,
};

const Args = struct {
    input_file: File,
    mines: MinesInfo,
    per_cell: u8 = 1,
    debug: bool = false,
    quiet: bool = false,
};

const CellContents = union(enum) {
    /// Unclicked cell.
    Unclicked,
    /// Revealed space (effectively number 0).
    Space,
    /// Revealed number (1, 2, 3, ...).
    Number: u8,
    /// Mine or flag (1, 2, 3, ...).
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

        pub fn init(x_size: usize, y_size: usize) !Self {
            const rows: [][]T = try allocator.alloc([]T, y_size);
            for (rows) |*row| {
                row.* = try allocator.alloc(T, x_size);
            }
            return Self{ .data = rows };
        }

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

        pub fn getEntryAtIdx(g: Self, idx: usize) Entry {
            const x = idx % g.xSize();
            const y = idx / g.xSize();
            return .{ .x = x, .y = y, .value = g.getCell(x, y) };
        }

        /// The 'idx' argument must be in range of the grid bounds, as given by
        /// 'Grid.ySize()'.
        pub fn getRow(g: Self, idx: usize) []const T {
            return g.data[idx];
        }

        /// The 'idx' argument must be in range of the grid bounds, as given by
        /// 'Grid.xSize()'. Allocates memory for the slice - memory owned by the
        /// caller.
        pub fn getColumn(g: Self, idx: usize) error{OutOfMemory}![]const T {
            var list = ArrayList(T).init(allocator);
            defer list.deinit();
            try list.ensureTotalCapacity(g.xSize());
            for (g.data) |row| {
                list.appendAssumeCapacity(row[idx]);
            }
            return list.toOwnedSlice();
        }

        /// The 'idx' argument must be in range of the grid bounds, as given by
        /// 'xSize()'. The provided slice must be big enough to store the column.
        pub fn getColumnIntoSlice(self: Self, idx: usize, slice: []T) void {
            const x = idx;
            var y: usize = 0;
            while (y < self.ySize()) : (y += 1) {
                slice[y] = self.getCell(x, y);
            }
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

        pub fn toStr(
            g: Self,
            comptime fmt: []const u8,
            opts: struct { sep_idx: ?usize = null },
        ) ![]const u8 {
            // Find the max cell width.
            var max_width: u64 = 0;
            var iter = g.iterator();
            while (iter.next()) |entry| {
                const cell_width = std.fmt.count(fmt, .{entry.value});
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
                    const cell_width = std.fmt.count(fmt, .{cell});
                    try writer.writeByteNTimes(' ', @intCast(usize, max_width - cell_width));
                    try writer.print(fmt, .{cell});
                }
            }
            return buf.toOwnedSlice();
        }

        /// Make a shallow copy of the grid - memory owned by the caller.
        pub fn copy(g: Self) error{OutOfMemory}!Self {
            return Self.fromOwnedData(g.data);
        }

        pub fn iterator(g: Self) Iterator {
            return .{ .grid = g };
        }

        /// Stores the result in the provided buffer.
        pub fn getNeighbours(g: Self, x: usize, y: usize, buf: *[8]Entry) []const Entry {
            var x_min: usize = x;
            var x_max: usize = x;
            var y_min: usize = y;
            var y_max: usize = y;
            if (x > 0) x_min -= 1;
            if (x < g.xSize() - 1) x_max += 1;
            if (y > 0) y_min -= 1;
            if (y < g.ySize() - 1) y_max += 1;

            var i: usize = x_min;
            var nbr_idx: usize = 0;
            while (i <= x_max) : (i += 1) {
                var j: usize = y_min;
                while (j <= y_max) : (j += 1) {
                    if (i == x and j == y) continue;
                    buf[nbr_idx] = .{ .x = i, .y = j, .value = g.getCell(i, j) };
                    nbr_idx += 1;
                }
            }
            return buf[0..nbr_idx];
        }

        /// Allocates the slice - memory owned by the caller.
        pub fn getNeighboursAlloc(g: Self, x: usize, y: usize) error{OutOfMemory}![]const Entry {
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
            while (i <= x_max) : (i += 1) {
                var j: usize = y_min;
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

    pub fn init(cells: [][]isize) error{OutOfMemory}!Self {
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
    pub fn selectColumns(self: Self, col_idxs: []const usize) error{OutOfMemory}!Self {
        const rows = try allocator.alloc([]isize, self.ySize());
        errdefer {
            for (rows) |row| allocator.free(row);
            allocator.free(rows);
        }
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
    /// 'xSize()'. Allocates memory for the slice - memory owned by the caller.
    pub fn getColumn(self: Self, idx: usize) error{OutOfMemory}![]const isize {
        return self.grid.getColumn(idx);
    }

    /// The 'idx' argument must be in range of the grid bounds, as given by
    /// 'xSize()'. The provided slice must be big enough to store the column.
    pub fn getColumnIntoSlice(self: Self, idx: usize, slice: []isize) void {
        self.grid.getColumnIntoSlice(idx, slice);
    }

    pub fn toStr(self: Self, opts: struct { sep_idx: ?usize = null }) ![]const u8 {
        return self.grid.toStr("{}", .{ .sep_idx = opts.sep_idx });
    }

    /// Returns a copy of the matrix - memory owned by the caller.
    pub fn copy(self: Self) error{OutOfMemory}!Self {
        return Self{ .grid = try self.grid.copy() };
    }

    pub fn matrixAdd(self: Self, other: Self) error{OutOfMemory}!Self {
        assert(self.xSize() == other.xSize() and self.ySize() == other.ySize());
        var new_matrix = try self.copy();
        for (new_matrix.grid.data) |row, y| {
            for (row) |*cell, x| {
                cell.* += other.getCell(x, y);
            }
        }
        return new_matrix;
    }

    pub fn matrixSubtract(self: Self, other: Self) error{OutOfMemory}!Self {
        assert(self.xSize() == other.xSize() and self.ySize() == other.ySize());
        var new_matrix = try self.copy();
        for (new_matrix.grid.data) |row, y| {
            for (row) |*cell, x| {
                cell.* -= other.getCell(x, y);
            }
        }
        return new_matrix;
    }

    pub fn matrixMultiply(self: Self, other: Self) error{OutOfMemory}!Self {
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

                const pivot_val = @intCast(usize, self.getCell(col_idx, row_idx));
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
            cell.* = @intCast(i16, cell.*) * @intCast(i16, factor);
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

    const Self = @This();

    /// Allocates memory, but also reuses the provided memory storing the
    /// 'cells' slice, which must all be managed by caller.
    pub fn fromFlatSlice(x_size: usize, y_size: usize, cells: []CellContents) !Self {
        return Self{
            .grid = try Grid(CellContents).fromFlatSlice(x_size, y_size, cells),
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
        return self.grid.toStr("{}", .{});
    }
};

const RectangularIterator = struct {
    max_vals: []const u16,
    /// Field to store memory to use for the return slice.
    result_slice: []u16,
    idx: usize = 0,

    pub fn init(max_vals: []const u16, result_slice: []u16) @This() {
        assert(result_slice.len == max_vals.len);
        return .{ .max_vals = max_vals, .result_slice = result_slice };
    }

    pub fn next(self: *@This()) ?[]const u16 {
        if (self.idx >= self.size()) return null;
        defer self.idx += 1;
        var counter: usize = 1;
        for (self.max_vals) |val, i| {
            self.result_slice[i] = @intCast(u16, (self.idx / counter) % (val + 1));
            counter *= val + 1;
        }
        return self.result_slice;
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
    mines: MinesInfo,
    per_cell: u8,
    computed_state: ComputedState = .{},

    const ComputedState = struct {
        /// Each inner slice contains indexes to grid positions.
        groups: ?[]const []const usize = null,
        /// Numbers displayed on the board, corresponding to rows in the matrix.
        numbers: ?[]const Number = null,
        /// Matrix of simultaneous equations corresponding to the board.
        matrix: ?Matrix = null,
        /// Each inner slice contains a number of mines for each group.
        configs: ?[]const []const u16 = null,

        const Number = struct {
            idx: usize,
            value: u8,
            effective_value: u8,
            groups: ArrayList(usize),
        };
    };

    const Self = @This();

    /// The board is still owned by the caller and should be independently
    /// deinitialised.
    pub fn init(board: Board, mines_info: MinesInfo, per_cell: u8) Self {
        return Self{
            .board = board,
            .mines = mines_info,
            .per_cell = per_cell,
        };
    }

    pub fn deinit(self: *Self) void {
        if (self.computed_state.groups) |groups| {
            for (groups) |grp| {
                allocator.free(grp);
            }
            allocator.free(groups);
            self.computed_state.groups = null;
        }

        if (self.computed_state.numbers) |numbers| {
            for (numbers) |num| {
                num.groups.deinit();
            }
            allocator.free(numbers);
            self.computed_state.numbers = null;
        }

        if (self.computed_state.matrix) |matrix| {
            matrix.deinit();
            self.computed_state.matrix = null;
        }

        if (self.computed_state.configs) |configs| {
            for (configs) |cfg| {
                allocator.free(cfg);
            }
            allocator.free(configs);
            self.computed_state.configs = null;
        }
    }

    /// High-level function to perform all steps required to solve the board,
    /// returning a grid of probabilities.
    pub fn solve(self: *Self) !Grid(f64) {
        if (self.computed_state.matrix == null or self.computed_state.groups == null) {
            try self.prepare();
            assert(self.computed_state.matrix != null);
            assert(self.computed_state.groups != null);
        }
        self.computed_state.configs = try self.findConfigs();
        return self.calcProbabilities();
    }

    /// High-level function to perform all preparation steps such as finding the
    /// matrix and equivalence groups.
    pub fn prepare(self: *Self) !void {
        std.log.info("Preparing solver with {d} x {d} board", .{ self.board.xSize(), self.board.ySize() });
        const cs = &self.computed_state;
        cs.matrix = try self.findGroupsMatrix();
        std.log.info("Initial matrix is {d} x {d}", .{ cs.matrix.?.xSize(), cs.matrix.?.ySize() });
        cs.matrix.?.rref();
        std.log.info("Reduced matrix to RREF", .{});
        cs.matrix.?.removeTrailingZeroRows();
        std.log.info("Removed zero-rows from matrix, {d} rows", .{cs.matrix.?.ySize()});
    }

    /// Convert a board into a set of simultaneous equations, represented in matrix
    /// form. The memory is owned by the caller.
    pub fn findFullMatrix(self: Self) !Matrix {
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
        const num_mines: ?u16 = if (self.mines == .Num) self.mines.Num else null;
        const num_columns = self.board.grid.count(.Unclicked) + 1;
        const num_rows = self.board.grid.count(.Number) + @as(usize, if (num_mines) |_| 1 else 0);
        const all_matrix_cells: []isize = try allocator.alloc(isize, num_columns * num_rows);
        defer allocator.free(all_matrix_cells);

        var j: usize = 0; // Matrix row index
        var iter1 = self.board.grid.iterator();
        // Iterate over the numbers in the board, which correspond to the rows of
        // the matrix.
        while (iter1.next()) |num_entry| {
            if (num_entry.value != .Number) continue;
            const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
            var i: usize = 0; // Matrix column index
            var iter2 = self.board.grid.iterator();
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

        // Add the final row on the end, corresponding to the total number of mines.
        if (num_mines) |m| {
            assert(j == num_rows - 1);
            var mines: u16 = 0;
            iter1.idx = 0;
            while (iter1.next()) |cell_entry| {
                if (cell_entry.value == .Mine) mines += cell_entry.value.Mine;
            }
            if (mines > m) return error.TooManyMinesInBoard;
            const row = all_matrix_cells[j * num_columns .. (j + 1) * num_columns];
            var i: usize = 0;
            while (i < num_columns - 1) : (i += 1) {
                row[i] = 1;
            }
            row[i] = m - mines;
        }

        return Matrix.fromFlatSlice(isize, num_columns, num_rows, all_matrix_cells);
    }

    /// Convert a board into a set of simultaneous equations for equivalence
    /// groups, represented in matrix form. The memory is owned by the caller.
    pub fn findGroupsMatrix(self: *Self) !Matrix {
        if (self.computed_state.groups == null) {
            self.computed_state.groups = try self.findGroups();
        }
        const groups = self.computed_state.groups.?;
        const numbers = self.computed_state.numbers.?;

        const num_columns = groups.len + 1;
        const num_rows = numbers.len + @as(usize, if (self.mines == .Num) 1 else 0);
        const all_cells = try allocator.alloc(isize, num_columns * num_rows);
        defer allocator.free(all_cells);
        for (all_cells) |*cell| cell.* = 0;
        const rows = try allocator.alloc([]isize, num_rows);
        defer allocator.free(rows);
        for (rows) |*row, j| row.* = all_cells[j * num_columns .. (j + 1) * num_columns];

        // Iterate over numbers, creating a matrix row for each.
        for (numbers) |num, j| {
            const row = rows[j];
            for (num.groups.items) |i| {
                row[i] = 1;
            }
            // Add the final column value - the RHS of the equation.
            row[row.len - 1] = num.effective_value;
        }
        // Add the final row on the end, corresponding to the total number of mines.
        if (self.mines == .Num) {
            const total_mines = self.mines.Num;
            assert(rows.len == numbers.len + 1);
            // Count all mines in the board.
            var disp_mines: u16 = 0;
            var iter1 = self.board.grid.iterator();
            while (iter1.next()) |cell_entry| {
                if (cell_entry.value == .Mine) disp_mines += cell_entry.value.Mine;
            }
            if (disp_mines > total_mines) return error.TooManyMinesInBoard;
            // Fill in the row.
            const row = all_cells[(rows.len - 1) * num_columns ..];
            var i: usize = 0;
            while (i < num_columns - 1) : (i += 1) {
                row[i] = 1;
            }
            row[i] = total_mines - disp_mines;
        }

        // TODO: Get contiguous block of memory working (double frees...)
        // return Matrix{ .grid = .{ .data = rows } };
        return Matrix.fromFlatSlice(isize, num_columns, num_rows, all_cells);
    }

    pub fn findGroups(self: *Self) ![]const []const usize {
        // - Iterate over all cells in the board
        // - Skip over cells that aren't unclicked (prob will be zero)
        // - Find the indices of neighbouring displayed numbers
        // - Check whether this combination of nbr nums has been found already
        //    - If it has, associate this cell with the corresponding group ID
        //    - Otherwise create a new group and associate this cell with it
        // Should end up with:
        //  - An array of groups (where elements correspond to board cell
        //    positions)
        //  - An array of numbers (where elements contain the cell index, the
        //    number value and a list of group indices)

        var groups = ArrayList(ArrayList(usize)).init(allocator);
        defer {
            for (groups.items) |g| g.deinit();
            groups.deinit();
        }
        var group_num_idxs = ArrayList([]const usize).init(allocator);
        defer {
            for (group_num_idxs.items) |g| allocator.free(g);
            group_num_idxs.deinit();
        }
        var numbers = ArrayList(ComputedState.Number).init(allocator);
        defer {
            for (numbers.items) |n| n.groups.deinit();
            numbers.deinit();
        }

        var nbrs_buf: [8]Grid(CellContents).Entry = undefined;
        var num_nbr_idxs_buf: [8]usize = undefined;
        var cell_idx: usize = 0;
        while (cell_idx < self.board.xSize() * self.board.ySize()) : (cell_idx += 1) {
            const entry = self.board.grid.getEntryAtIdx(cell_idx);
            if (entry.value != .Unclicked) continue;

            // Create a slice of indices of number neighbours.
            const nbrs = self.board.grid.getNeighbours(entry.x, entry.y, &nbrs_buf);
            var num_nbr_idxs_idx: usize = 0;
            for (nbrs) |nbr| {
                if (nbr.value == .Number) {
                    num_nbr_idxs_buf[num_nbr_idxs_idx] = nbr.x + self.board.xSize() * nbr.y;
                    num_nbr_idxs_idx += 1;
                }
            }
            const num_nbr_idxs = num_nbr_idxs_buf[0..num_nbr_idxs_idx];

            // Don't include outer group if using mine density.
            if (num_nbr_idxs.len == 0 and self.mines == .Density) continue;

            // Check whether this combination of numbers has already been found.
            var grp_idx: usize = 0;
            while (grp_idx < groups.items.len) : (grp_idx += 1) {
                if (std.mem.eql(usize, num_nbr_idxs, group_num_idxs.items[grp_idx])) {
                    break;
                }
            } else {
                // Start a new group.
                try groups.append(ArrayList(usize).init(allocator));
                try group_num_idxs.append(try allocator.dupe(usize, num_nbr_idxs));
                for (num_nbr_idxs) |num_idx| {
                    // Check whether this number has already been found.
                    var i: usize = 0;
                    while (i < numbers.items.len) : (i += 1) {
                        if (numbers.items[i].idx == num_idx) {
                            break;
                        }
                    } else {
                        const num_entry = self.board.grid.getEntryAtIdx(num_idx);
                        const value = num_entry.value.Number;
                        // Count the neighbouring mines.
                        var mines: u16 = 0;
                        // Note: reusing nbrs_buf - be careful that this is valid!
                        const num_nbrs = self.board.grid.getNeighbours(
                            num_entry.x,
                            num_entry.y,
                            &nbrs_buf,
                        );
                        for (num_nbrs) |num_nbr_entry| {
                            if (num_nbr_entry.value == .Mine)
                                mines += num_nbr_entry.value.Mine;
                        }
                        if (mines > value) return error.TooManyMinesAroundNumber;
                        try numbers.append(.{
                            .idx = num_idx,
                            .value = value,
                            .effective_value = value - @intCast(u8, mines),
                            .groups = ArrayList(usize).init(allocator),
                        });
                    }
                    try numbers.items[i].groups.append(grp_idx);
                }
            }
            try groups.items[grp_idx].append(cell_idx);
        }

        self.computed_state.numbers = numbers.toOwnedSlice();

        const groups_slice = try allocator.alloc([]const usize, groups.items.len);
        errdefer allocator.free(groups_slice);
        for (groups.items) |*grp, i| {
            groups_slice[i] = grp.toOwnedSlice();
        }
        return groups_slice;
    }

    /// Returns a slice containing valid mine configurations, which are
    /// represented as slices of length equal to the number of groups (i.e. the
    /// number of matrix columns), with the values being the number of mines in
    /// the corresponding group.
    ///
    /// Memory owned by the caller.
    fn findConfigs(self: Self) ![]const []const u16 {
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
        // These equations are subjust to g_i >= 0, which allows us to constrain
        // the search space with e.g.:
        //   |  0  1  1 |             |  1 |
        //   |  1  0  0 |   | g6 |    |  4 |
        //   |  0 -1 -1 | . | g7 | <= |  1 |
        //   | -1  1  0 |   | g8 |    | -2 |
        //   |  1  0  1 |             |  4 |
        //

        const matrix = self.computed_state.matrix.?;
        const groups = self.computed_state.groups.?;

        // Check for inconsistent matrix - row with all-zero LHS, non-zero RHS
        // (last row in RREF with zero-rows removed).
        {
            const last_row = matrix.grid.data[matrix.ySize() - 1];
            if (std.mem.allEqual(isize, last_row[0 .. matrix.xSize() - 1], 0) and
                last_row[matrix.xSize() - 1] != 0)
            {
                return error.InvalidMatrixEquations;
            }
        }

        const rhs_vec = try matrix.selectColumns(&[_]usize{matrix.xSize() - 1});
        defer rhs_vec.deinit();
        const col_categorisation = try self.categoriseColumns();
        const fixed_col_idxs = col_categorisation.fixed;
        defer allocator.free(fixed_col_idxs);
        const free_col_idxs = col_categorisation.free;
        defer allocator.free(free_col_idxs);
        std.log.info(
            "Categorised columns: {d} fixed, {d} free",
            .{ fixed_col_idxs.len, free_col_idxs.len },
        );

        var configs = ArrayList([]u16).init(allocator);
        errdefer for (configs.items) |cfg| allocator.free(cfg);
        defer configs.deinit();

        if (free_col_idxs.len == 0) {
            const config = try allocator.alloc(u16, groups.len);
            errdefer allocator.free(config);
            for (config) |*val, i| {
                if (rhs_vec.getCell(0, i) < 0) return error.InvalidMatrixEquations;
                val.* = @intCast(u16, rhs_vec.getCell(0, i));
            }
            std.log.info("Found a unique solution", .{});
            try configs.append(config);
            return configs.toOwnedSlice();
        }

        const free_var_matrix = try matrix.selectColumns(free_col_idxs);

        const max_free_vals = try self.findMaxFreeVals(free_var_matrix, rhs_vec, free_col_idxs);
        defer allocator.free(max_free_vals);
        const free_vals_storage = try allocator.alloc(u16, max_free_vals.len);
        defer allocator.free(free_vals_storage);
        var iter = RectangularIterator.init(max_free_vals, free_vals_storage);
        std.log.info(
            "Using free val maximums: {d} ({d} combinations)",
            .{ max_free_vals, iter.size() },
        );

        const fixed_vals = try allocator.alloc(isize, fixed_col_idxs.len);
        defer allocator.free(fixed_vals);
        while (iter.next()) |free_vals| {
            // Fill in 'fixed_vals' by performing matrix multiplication.
            for (fixed_vals) |*cell, i| {
                var subtract_val: isize = 0;
                for (free_var_matrix.getRow(i)) |free_cell, j| {
                    subtract_val += free_cell * free_vals[j];
                }
                cell.* = rhs_vec.getCell(0, i) - subtract_val;
            }

            var invalid_var_idx: ?usize = null;
            for (fixed_vals) |cell, idx| {
                if (cell < 0 or
                    cell > groups[fixed_col_idxs[idx]].len * self.per_cell)
                {
                    invalid_var_idx = idx;
                    break;
                }
            }
            if (invalid_var_idx) |idx| {
                std.log.debug(
                    "Potential config {d} invalid in fixed var {d}",
                    .{ iter.idx, idx },
                );
                continue;
            }

            var config = try allocator.alloc(u16, groups.len);
            errdefer allocator.free(config);

            for (fixed_col_idxs) |grp_idx, fixed_col_idx| {
                config[grp_idx] = @intCast(u16, fixed_vals[fixed_col_idx]);
            }
            for (free_col_idxs) |grp_idx, free_col_idx| {
                config[grp_idx] = @intCast(u16, free_vals[free_col_idx]);
            }
            try configs.append(config);
            std.log.debug(
                "Potential config {d} is valid number {d}",
                .{ iter.idx, configs.items.len },
            );
        }

        std.log.info("Found {d} mine configurations", .{configs.items.len});

        if (configs.items.len == 0) return error.InvalidMatrixEquations;

        return configs.toOwnedSlice();
    }

    fn calcProbabilities(self: Self) !Grid(f64) {
        const groups = self.computed_state.groups.?;
        const configs = self.computed_state.configs.?;

        std.log.debug("Calculating config probs...", .{});
        const cfg_probs = try allocator.alloc(f64, configs.len);
        for (configs) |cfg, idx| {
            var log_combs: f64 = 0;
            for (cfg) |m_i, i| {
                const g_size = groups[i].len;
                log_combs += logCombs(g_size, m_i, self.per_cell);
                var k: u16 = 1;
                while (k <= m_i) : (k += 1) { // Divide by m_i!
                    log_combs -= std.math.ln(@intToFloat(f64, k));
                }
                switch (self.mines) { // Multiply by  ( rho / (1-rho) )^m_i
                    .Density => |rho| {
                        log_combs += std.math.ln(rho / (1 - rho)) * @intToFloat(f64, m_i);
                    },
                    else => {},
                }
            }
            cfg_probs[idx] = std.math.exp(log_combs);
            assert(cfg_probs[idx] > 0);
        }

        var weight: f64 = 0;
        for (cfg_probs) |p| weight += p;
        if (weight == std.math.inf(f64)) return error.TooManyCombinations;
        std.log.debug("Adjusting config probs by weight {e:.2}", .{weight});
        for (cfg_probs) |*p, i| {
            p.* = p.* / weight;
            std.log.debug("Prob for config {d}: {d:.4}", .{ i, p.* });
            assert(p.* <= 1.0001);
        }

        const x_size = self.board.xSize();
        const y_size = self.board.ySize();
        const probs_grid = try Grid(f64).init(x_size, y_size);
        for (probs_grid.data) |row| {
            for (row) |*cell| {
                cell.* = 0;
            }
        }

        std.log.debug("Calculating group probs from config probs...", .{});
        for (groups) |grp_i, i| {
            var unsafe_prob: f64 = 0;
            for (configs) |cfg_j, j| {
                unsafe_prob += cfg_probs[j] * unsafeProb(grp_i.len, cfg_j[i], self.per_cell);
            }
            std.log.debug("Calculated unsafe prob for group {d}: {d:.3}", .{ i, unsafe_prob });
            assert(unsafe_prob >= 0);
            // TODO: Add in some rounding to ensure no higher than 1!
            assert(unsafe_prob <= 1.0001);
            for (grp_i) |cell_idx| {
                const entry = probs_grid.getEntryAtIdx(cell_idx);
                probs_grid.data[entry.y][entry.x] = unsafe_prob;
            }
        }

        // TODO: If working with density then we need to fill in the outer group.

        return probs_grid;
    }

    const ColumnCategorisation = struct { fixed: []const usize, free: []const usize };

    fn categoriseColumns(self: Self) !ColumnCategorisation {
        const matrix = self.computed_state.matrix.?;

        var fixed_cols = ArrayList(usize).init(allocator);
        defer fixed_cols.deinit();
        var free_cols = ArrayList(usize).init(allocator);
        defer free_cols.deinit();

        var x: usize = 0;
        var y: usize = 0;
        while (y < matrix.ySize()) : (y += 1) {
            while (x < matrix.xSize() - 1) : (x += 1) {
                if (matrix.getCell(x, y) == 0) {
                    try free_cols.append(x);
                } else {
                    try fixed_cols.append(x);
                    x += 1;
                    break;
                }
            }
        }
        // All remaining columns are free columns (leading '1's all found).
        while (x < matrix.xSize() - 1) : (x += 1) {
            try free_cols.append(x);
        }
        return ColumnCategorisation{
            .fixed = fixed_cols.toOwnedSlice(),
            .free = free_cols.toOwnedSlice(),
        };
    }

    fn findMaxFreeVals(
        self: Self,
        free_var_matrix: Matrix,
        rhs_vec: Matrix,
        free_col_idxs: []const usize,
    ) error{OutOfMemory}![]const u16 {
        assert(free_var_matrix.ySize() == rhs_vec.ySize());
        assert(free_var_matrix.xSize() == free_col_idxs.len);
        assert(rhs_vec.xSize() == 1);

        const groups = self.computed_state.groups.?;

        const max_free_vals = try allocator.alloc(u16, free_var_matrix.xSize());
        for (max_free_vals) |*val, i| {
            val.* = @intCast(u16, groups[free_col_idxs[i]].len) * self.per_cell;
            var y: usize = 0;
            while (y < free_var_matrix.ySize()) : (y += 1) {
                if (rhs_vec.getCell(0, y) < 0) continue;
                if (free_var_matrix.getCell(i, y) <= 0) continue;
                var unusable_row = false;
                for (free_var_matrix.getRow(y)) |cell| {
                    if (cell < 0) unusable_row = true;
                }
                if (unusable_row) continue;

                const bound =
                    @intCast(u16, rhs_vec.getCell(0, y)) /
                    @intCast(u16, free_var_matrix.getCell(i, y));
                if (bound < val.*) val.* = bound;
            }
        }
        return max_free_vals;
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

fn parseInputBoard(input: []const u8) !Board {
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
    );
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

fn parseArgs() !Args {
    // First we specify what parameters our program can take.
    // We can use 'parseParam()' to parse a string to a 'Param(Help)'.
    const params = comptime [_]clap.Param(clap.Help){
        clap.parseParam("-h, --help                    Display this help and exit") catch unreachable,
        clap.parseParam("-f, --file <PATH>             Input file (defaults to stdin)") catch unreachable,
        clap.parseParam("-m, --mines <MINES>           Number of mines") catch unreachable,
        clap.parseParam("-d, --infinite-density <VAL>  Density of mines on infinite board") catch unreachable,
        clap.parseParam("-p, --per-cell <NUM>          Max number of mines per cell") catch unreachable,
        clap.parseParam("-v, --verbose                 Output debug info and logging to stderr") catch unreachable,
        clap.parseParam("-q, --quiet                   Emit less output to stderr") catch unreachable,
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

    const input_file = blk: {
        if (clap_args.option("--file")) |file| {
            break :blk std.fs.cwd().openFile(file, .{}) catch |err| {
                try stderr.print("Failed to open file {s}\n", .{file});
                return err;
            };
        } else break :blk std.io.getStdIn();
    };

    const mines_arg_set = clap_args.option("--mines") != null;
    const density_arg_set = clap_args.option("--infinite-density") != null;
    if (mines_arg_set and density_arg_set or !mines_arg_set and !density_arg_set) {
        try stderr.writeAll("Exactly one of '--mines' and '--infinite-density' args expected\n");
        return error.InvalidArgument;
    }

    const mines: MinesInfo = blk: {
        if (clap_args.option("--mines")) |num_mines_str| {
            const num_mines = std.fmt.parseUnsigned(u16, num_mines_str, 10) catch |err| {
                try stderr.writeAll("Expected positive integer number of mines\n");
                return err;
            };
            if (num_mines == 0) {
                try stderr.writeAll("Expected positive integer number of mines\n");
                return error.InvalidArgument;
            }
            break :blk .{ .Num = num_mines };
        } else if (clap_args.option("--infinite-density")) |density_str| {
            const density = std.fmt.parseFloat(f64, density_str) catch |err| {
                try stderr.writeAll("Expected mines density between 0 and 1\n");
                return err;
            };
            if (density <= 0 or density >= 1) {
                try stderr.writeAll("Expected mines density between 0 and 1\n");
                return error.InvalidArgument;
            }
            break :blk .{ .Density = density };
        } else unreachable;
    };

    var args = Args{ .input_file = input_file, .mines = mines };

    if (clap_args.option("--per-cell")) |per_cell| {
        args.per_cell = try std.fmt.parseUnsigned(u8, per_cell, 10);
        if (args.per_cell == 0) {
            try stderr.writeAll("Max number of mines per cell must be greater than 0\n");
            return error.InvalidArgument;
        }
    }

    if (clap_args.flag("--verbose")) {
        args.debug = true;
    }
    if (clap_args.flag("--quiet")) {
        args.quiet = true;
    }
    if (args.quiet and args.debug) {
        try stderr.writeAll("--quiet and --verbose cannot be specified together\n");
        return error.InvalidArgument;
    }

    return args;
}

pub fn main() !u8 {
    start_milli_time = @intCast(u64, std.time.milliTimestamp());
    // Set up an allocator - no need to free memory as we go.
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    allocator = &gpa.allocator;

    stdout = std.io.getStdOut().writer();
    stderr = std.io.getStdErr().writer();

    const args = parseArgs() catch |err| switch (err) {
        error.InvalidArgument,
        error.MissingValue,
        error.DoesntTakeValue,
        => return 2,
        else => return 1,
    };
    if (args.debug) {
        internal_log_level = .debug;
    } else if (args.quiet) {
        internal_log_level = .notice;
    }

    // Debug log the parsed args.
    const verbosity = if (args.debug) "verbose" else if (args.quiet) "quiet" else "default";
    switch (args.mines) {
        .Num => |num_mines| std.log.debug(
            "Parsed args: mines={d}, per_cell={d}, verbosity={s}",
            .{ num_mines, args.per_cell, verbosity },
        ),
        .Density => |density| std.log.debug(
            "Parsed args: density={d}, per_cell={d}, verbosity={s}",
            .{ density, args.per_cell, verbosity },
        ),
    }

    if (args.per_cell > 3) {
        std.log.err(
            "Currently unable to perform probability calculation for per_cell > 3",
            .{},
        );
        return 1;
    }

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

    const board = try parseInputBoard(input);
    defer board.deinit();
    try stdout.print("Board:\n{s}\n", .{try board.toStr()});

    var solver = Solver.init(board, args.mines, args.per_cell);
    defer solver.deinit();

    if (args.debug) {
        // var matrix = try solver.findFullMatrix();
        // defer matrix.deinit();
        // std.log.info("", .{});
        // std.log.info(
        //     "Matrix:\n{s}",
        //     .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
        // );

        // matrix.rref();
        // std.log.info("", .{});
        // std.log.info(
        //     "RREF matrix:\n{s}",
        //     .{try matrix.toStr(.{ .sep_idx = matrix.xSize() - 1 })},
        // );

        const groups = try solver.findGroups();
        const cs = &solver.computed_state;

        std.log.info("", .{});
        std.log.info("Solver numbers:", .{});
        for (cs.numbers.?) |num, i| {
            std.log.info("{d}: {{.idx = {d}, .value = {d}}}", .{ i, num.idx, num.value });
        }

        std.log.info("", .{});
        std.log.info("Solver groups:", .{});
        for (groups) |grp, i| {
            std.log.info("{d}: {d}", .{ i, grp });
        }

        const grp_matrix = try solver.findGroupsMatrix();

        std.log.info("", .{});
        std.log.info(
            "Solver groups matrix:\n{s}",
            .{try grp_matrix.toStr(.{ .sep_idx = grp_matrix.xSize() - 1 })},
        );

        try solver.prepare();

        std.log.info("", .{});
        std.log.info(
            "Solver final matrix:\n{s}",
            .{try cs.matrix.?.toStr(.{ .sep_idx = cs.matrix.?.xSize() - 1 })},
        );
    }

    const probs = try solver.solve();

    if (args.debug) {
        std.log.info("", .{});
        std.log.info("Mine configurations:", .{});
        for (solver.computed_state.configs.?) |cfg, i| {
            std.log.info("{d}: {d}", .{ i, cfg });
        }
    }

    try stdout.print("\n", .{});

    try stdout.print("Probabilities:\n{s}\n", .{try probs.toStr("{d:.5}", .{})});

    std.log.info("Finished", .{});
    return 0;
}

// -----------------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------------

test "lcm" {
    try std.testing.expectEqual(@as(usize, 2), lcm(1, 2));
    try std.testing.expectEqual(@as(usize, 10), lcm(5, 2));
    try std.testing.expectEqual(@as(usize, 24), lcm(8, 6));
}

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
    const max_vals = [_]u16{ 2, 0, 1 };
    var slice_mem = [_]u16{0} ** 3;
    var iter = RectangularIterator{ .max_vals = &max_vals, .result_slice = &slice_mem };
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
        if (iter.next()) |vals| {
            try std.testing.expectEqualSlices(u16, exp_vals, vals);
        } else return error.ExpectedIteratorItem;
    }
    try std.testing.expectEqual(@as(?[]const u16, null), iter.next());
}

test "Solver: invalid board" {
    allocator = std.testing.allocator;
    const board = try parseInputBoard(
        \\ # 1 2
        \\ # # #
    );
    defer board.deinit();

    var solver = Solver.init(board, .{ .Num = 1 }, 1);
    defer solver.deinit();

    try std.testing.expectError(error.InvalidMatrixEquations, solver.solve());
}
