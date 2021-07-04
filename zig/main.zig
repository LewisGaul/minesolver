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
};

const CellContents = union(enum) {
    Unclicked,
    Num: u8,
    Mine: u8,

    pub fn format(
        self: @This(),
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        switch (self) {
            .Unclicked => try writer.writeByte('#'),
            .Num => |n| {
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
        x_size: u8,
        y_size: u8,
        cells: []T,

        const Self = @This();

        /// Memory owned by callee, i.e. ensure the slice of cells will live
        /// for long enough!
        pub fn init(x_size: u8, y_size: u8, cells: []T) !Self {
            if (cells.len != x_size * y_size) return error.InvalidNumberOfCells;
            return Self{ .x_size = x_size, .y_size = y_size, .cells = cells };
        }

        pub fn get(g: Self, x: u8, y: u8) !T {
            if (x >= g.x_size or y >= g.y_size) return error.OutOfBounds;
            return g.cells[y * g.x_size + x];
        }

        pub fn toStr(g: Self) ![]const u8 {
            var buf = ArrayList(u8).init(allocator);
            errdefer buf.deinit();
            var writer = buf.writer();
            for (g.cells) |c, i| {
                if (i > 0 and i % g.x_size == 0) try writer.writeByte('\n');
                try writer.print("{} ", .{c});
            }
            return buf.toOwnedSlice();
        }
    };
}

const Board = Grid(CellContents);

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

fn parseInputCell(input: []const u8) !CellContents {
    switch (input[0]) {
        '0'...'9' => return CellContents{
            .Num = try std.fmt.parseUnsigned(u8, input, 10),
        },
        '.' => {
            if (input.len > 1) return error.UnexpectedCellText;
            return CellContents{ .Num = 0 };
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

    return Board.init(first_line_cols.?, rows, cells_array.toOwnedSlice());
}

fn parseArgs() !Args {
    // First we specify what parameters our program can take.
    // We can use 'parseParam()' to parse a string to a 'Param(Help)'.
    const params = comptime [_]clap.Param(clap.Help){
        clap.parseParam("-h, --help           Display this help and exit") catch unreachable,
        clap.parseParam("-f, --file <PATH>    Input file (defaults to stdin)") catch unreachable,
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

    return Args{ .input_file = input_file };
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

    return 0;
}
