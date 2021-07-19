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

// -----------------------------------------------------------------------------
// Declarations
// -----------------------------------------------------------------------------

var allocator: *Allocator = undefined;

const stderr = std.io.getStdErr().writer();

const Args = struct {
    input_file: File,
    mines: u8,
    per_cell: u8 = 1,
};

const CellContents = union(enum) {
    Unclicked,
    Number: u8,
    Mine: u8,

    pub fn format(
        self: @This(),
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        switch (self) {
            .Unclicked => try writer.writeByte('#'),
            .Number => |n| {
                if (n == 0) {
                    try writer.writeByte('.');
                } else {
                    try writer.print("{}", .{n});
                }
            },
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

        /// Iterator for iterating over cells in index order.
        const Iterator = struct {
            grid: Grid,
            idx: u8 = 0,

            const Entry = struct { x: u8, y: u8, value: T };

            pub fn next(it: *@This()) ?Entry {
                if (it.idx >= it.grid.xSize() * it.grid.ySize()) return null;
                const x = it.idx % it.grid.ySize();
                const y = @divFloor(it.idx, it.grid.ySize());
                const entry = Entry{
                    .x = x,
                    .y = y,
                    .value = it.grid.get(x, y) catch unreachable,
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
            return g.data[0].len;
        }

        pub fn ySize(g: Self) u8 {
            return g.data.len;
        }

        pub fn get(g: Self, x: u8, y: u8) !T {
            if (x >= g.xSize() or y >= g.ySize()) return error.OutOfBounds;
            return g.data[y][x];
        }

        pub fn count(g: Self, item: T) usize {
            var total: usize = 0;
            for (g.data) |row| {
                total += std.mem.count(T, row, &.{item});
            }
            return total;
        }

        pub fn toStr(g: Self) ![]const u8 {
            var buf = ArrayList(u8).init(allocator);
            errdefer buf.deinit();
            var writer = buf.writer();
            for (g.data) |row, y| {
                if (y > 0) try writer.writeByte('\n');
                for (row) |cell| {
                    try writer.print("{} ", .{cell});
                }
            }
            return buf.toOwnedSlice();
        }

        pub fn iterator(g: Self) Iterator {
            return .{ .grid = g };
        }
    };
}

const Board = Grid(CellContents);
const Matrix = Grid(u8);

// -----------------------------------------------------------------------------
// Matrix board representation
// -----------------------------------------------------------------------------

/// Convert a board into a set of simultaneous equations, represented in matrix
/// form. The memory is owned by the caller.
fn boardToMatrix(board: Board) !Matrix {
    const unclicked_cells = 10;
    // const unclicked_cells = board.count(CellContents.Unclicked); TODO
    const number_cells = 5;
    // const number_cells = board.count(CellContents.Number); TODO
    const rows: [][]u8 = try allocator.alloc([]u8, number_cells);
    const all_matrix_cells: []u8 = try allocator.alloc(u8, unclicked_cells * number_cells);
    var y: u8 = 0;
    while (y < number_cells) : (y += 1) {
        rows[y] = all_matrix_cells[y * unclicked_cells .. (y + 1) * unclicked_cells];
    }
    return Matrix{ .data = rows };
}

// -----------------------------------------------------------------------------
// Input parsing
// -----------------------------------------------------------------------------

fn parseInputCell(input: []const u8) !CellContents {
    switch (input[0]) {
        '0'...'9' => return CellContents{
            .Number = try std.fmt.parseUnsigned(u8, input, 10),
        },
        '.' => {
            if (input.len > 1) return error.UnexpectedCellText;
            return CellContents{ .Number = 0 };
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

    return Board.fromFlatSlice(first_line_cols.?, rows, cells_array.toOwnedSlice());
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
    const mines = std.fmt.parseUnsigned(u8, mines_arg, 10) catch |err| {
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

    const board = try parseInputBoard(input);
    std.debug.print("{s}\n", .{try board.toStr()});

    const matrix = try boardToMatrix(board);
    std.debug.print("{s}\n", .{try matrix.toStr()});

    return 0;
}
