//
// Take an input minesweeper board and output info on the solving calculations.
//
// Input format:
//  - Each non-blank line corresponds to a row
//  - Cells in a row are separated by any number of spaces or tabs
//  - A hash symbol (#) represents an unclicked cell
//  - Numbers represent numbers shown on the board
//  - An asterisk (*) represents a mine, and may optionally be followed by a
//    number to indicate the number of mines (1 assumed otherwise)
//

const std = @import("std");

const clap = @import("clap");

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

        pub fn init(gpa: *Allocator, x_size: u8, y_size: u8, cells: []T) !Self {
            if (cells.len != x_size * y_size) return error.InvalidNumberOfCells;
            return Self{
                .x_size = x_size,
                .y_size = y_size,
                .cells = try gpa.dupe(T, cells),
            };
        }

        pub fn get(g: Self, x: u8, y: u8) !T {
            if (x >= g.x_size or y >= g.y_size) return error.OutOfBounds;
            return g.cells[y * g.x_size + x];
        }
    };
}

const Board = Grid(CellContents);

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------

fn parseInputBoard(input: []const u8) !Board {
    return Board.init(allocator, 0, 0, &.{});
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

    return 0;
}
