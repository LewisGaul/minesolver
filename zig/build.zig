const Builder = @import("std").build.Builder;

pub fn build(b: *Builder) void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Standard release options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall.
    const mode = b.standardReleaseOptions();

    // Building the main exe.
    const exe = b.addExecutable("zig-main", "main.zig");
    exe.setTarget(target);
    exe.setBuildMode(mode);
    exe.addPackagePath("clap", "deps/zig-clap/clap.zig");
    exe.install();

    // Running tests.
    var tests = b.addTest("main.zig");

    // Define the 'test' subcommand.
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&tests.step);

    const run_cmd = exe.run();
    run_cmd.step.dependOn(&exe.step); // TODO: Is this needed?

    // Define the 'run' subcommand.
    const run_step = b.step("run", "Run the script");
    run_step.dependOn(&run_cmd.step);
}
