import textwrap

import setuptools


long_description = textwrap.dedent(
    """\
    # Minesweeper solver

    A minesweeper solver written in Zig.


    The following cell representations are used for input boards:
    - `#` for an unclicked cell
    - `<N>` where `N=0,1,2,...` is a number shown in a cell
    - `.` as an alternative to `0` (since the number 0 is not normally shown)
    - `*` to represent a single mine (may be a revealed mine or a flag)
    - `*<N>` where `N=1,2,...` is the number of mines


    Example usage:
    ```
    >>> import zig_minesolver
    >>> board = \"""
    ... # 2 # # #
    ... # # # # #
    ... # 3 # # #
    ... # 2 # 4 #
    ... # # # # #
    ... \"""
    >>> probs = zig_minesolver.get_board_probs(board, mines=8)
    >>> print("\n".join(str(x) for x in probs))
    [0.27108, 0.0, 0.27108, 0.31325, 0.31325]
    [0.48594, 0.48594, 0.48594, 0.31325, 0.31325]
    [0.26506, 0.0, 0.50602, 0.5494, 0.5494]
    [0.26506, 0.0, 0.50602, 0.0, 0.5494]
    [0.10843, 0.10843, 0.24096, 0.5494, 0.5494]
    >>>
    ```
    """
)


# TODO: Currently unable to build as part of a pip install, reliant on
#  availability of wheels.


# class my_build(distutils.command.build.build):
#     """
#     Add custom steps into the build, which is executed when running 'pip install'.
#     """
#
#     def run(self):
#         compile_zig_project()
#         super().run()


setuptools.setup(
    name="zig-minesolver",
    version="0.1.1",
    author="Lewis Gaul",
    maintainer="Lewis Gaul",
    description="Minesweeper solver written in Zig",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LewisGaul/minesolver",
    # cmdclass={"build": my_build},
    packages=["minesolver"],
    package_data={"minesolver": ["zig-main*"]},
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Development Status :: 3 - Alpha",
    ],
    zip_safe=False,
)
