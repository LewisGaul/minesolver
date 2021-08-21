__all__ = ("Board", "CellContents", "Coord", "Grid")

import copy
from typing import Any, Callable, Collection, Generic, Iterable, Tuple, Type, TypeVar

import adt


Coord = Tuple[int, int]


class CellContents(metaclass=adt.ADTMeta):

    Unclicked: ()
    Num: (int,)
    Flag: (int,)
    WrongFlag: (int,)
    Mine: (int,)
    HitMine: (int,)

    @classmethod
    def from_string(cls, string: str) -> "CellContents":
        """
        Get the ADT field using the string representation.

        :param string:
            The string representation of a cell contents.
        :return:
            The cell contents.
        """
        try:
            if string.isnumeric():
                return cls.Num(int(string))
            elif string == ".":
                return cls.Num(0)
            elif string == "#":
                return cls.Unclicked()
            else:
                if len(string) == 2:
                    char, num = string
                    num = int(num)
                else:
                    assert len(string) == 1
                    char, num = string, 1
                if char == "F":
                    return cls.Flag(num)
                elif char == "X":
                    return cls.WrongFlag(num)
                elif char == "M":
                    return cls.Mine(num)
                elif char == "!":
                    return cls.HitMine(num)
                else:
                    assert False
        except Exception:
            raise ValueError(
                f"Unrecognised cell contents representation {string!r}"
            ) from None

    @adt.fieldmethod
    def to_string(field, basecls: Type["CellContents"]) -> str:
        """
        Convert an ADT field to the string representation.

        :return:
            The string representation.
        """
        if type(field) is basecls.Unclicked:
            return "#"
        elif type(field) is basecls.Num:
            if field[0] == 0:
                return "."
            else:
                return str(field[0])
        elif type(field) is basecls.Flag:
            return f"F{field[0]}"
        elif type(field) is basecls.WrongFlag:
            return f"X{field[0]}"
        elif type(field) is basecls.Mine:
            return f"M{field[0]}"
        elif type(field) is basecls.HitMine:
            return f"!{field[0]}"

    @adt.fieldmethod
    def is_mine_type(field, basecls: Type["CellContents"]) -> bool:
        """Determine whether a field is a mine type."""
        return field in [basecls.Flag, basecls.WrongFlag, basecls.Mine, basecls.HitMine]


T = TypeVar("T")


class Grid(Generic[T]):
    """
    Grid representation (square 2D array).

    Cells accessed by coordinate indexing, e.g. grid[x, y], where (0, 0) is
    the top-left corner.
    """

    def __init__(self, x_size: int, y_size: int, *, fill: T = 0):
        """
        :param x_size:
            The number of columns.
        :param y_size:
            The number of rows.
        :param fill:
            What to fill the grid with.
        """
        self._data = [x_size * [fill] for _ in range(y_size)]
        self.x_size: int = x_size
        self.y_size: int = y_size

    def __repr__(self):
        return f"<{self.x_size} x {self.y_size} {type(self).__name__}>"

    def __str__(self):
        # Use max length of object representation.
        cell_size = max(len(repr(obj)) for obj in self.values)
        cell_size = min(cell_size, 10)

        cell = f"{{:>{cell_size}}}"
        lines = []
        for row in self.rows:
            reprs = (cell.format(repr(obj)[:cell_size]) for obj in row)
            lines.append(" ".join(reprs))
        return "\n".join(lines)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if (self.x_size, self.y_size) != (other.x_size, other.y_size):
            return False
        for coord in self.coords:
            if self[coord] != other[coord]:
                return False
        return True

    def __contains__(self, coord: Coord):
        """Whether the given coord is in range of the grid."""
        try:
            x, y = coord
            return 0 <= x < self.x_size and 0 <= y < self.y_size
        except Exception:
            return False

    def __getitem__(self, key: Coord) -> T:
        try:
            x, y = key
            return self._data[y][x]
        except Exception:
            raise TypeError("Grid access expects a tuple coordinate of the form (0, 1)")

    def __setitem__(self, key: Coord, value: T):
        try:
            x, y = key
            self._data[y][x] = value
        except Exception:
            raise TypeError("Grid access expects a tuple coordinate of the form (0, 1)")

    @classmethod
    def from_grid(cls, grid: "Grid") -> "Grid":
        ret = cls(grid.x_size, grid.y_size)
        for coord in grid.coords:
            ret[coord] = grid[coord]
        return ret

    @classmethod
    def from_flat_iter(cls, iterable: Iterable[T], x_size: int, y_size: int) -> "Grid":
        """
        Create an instance from a flat iterable.

        :param iterable:
            The iterable to create the grid instance from. Must have a length
            matching the given dimensions.
        :param x_size:
            The number of columns.
        :param y_size:
            The number of rows.
        :return:
            The created grid.
        """
        grid = cls(x_size, y_size)
        i = 0
        for i, val in enumerate(iterable):
            x, y = i % x_size, i // x_size
            grid[x, y] = val
        if i != x_size * y_size - 1:
            raise ValueError(
                f"Size of iterable does not match dimensions: {i+1} != {x_size} x {y_size}"
            )
        return grid

    @classmethod
    def from_2d_iter(cls, iter_2d: Iterable[Iterable[T]]) -> "Grid":
        """
        Create an instance from a 2-dimensional array.

        Arguments:
        array ([[object, ...], ...])
            The array to use in creating the grid instance.

        Return: Grid
            The resulting grid.
        """
        try:
            x_size = len(iter_2d[0])
            y_size = len(iter_2d)
        except TypeError:
            iter_2d = [list(row) for row in iter_2d]
            x_size = len(iter_2d[0])
            y_size = len(iter_2d)
        return cls.from_flat_iter(
            (val for row in iter_2d for val in row), x_size, y_size
        )

    @property
    def coords(self) -> Iterable[Coord]:
        return ((x, y)  for y in range(self.y_size) for x in range(self.x_size))

    @property
    def values(self) -> Iterable[T]:
        return (v for row in self._data for v in row)

    @property
    def rows(self) -> Iterable[Iterable[T]]:
        return (iter(row) for row in self._data)

    @property
    def columns(self) -> Iterable[Iterable[T]]:
        return (
            (self._data[y][x] for y in range(self.y_size)) for x in range(self.x_size)
        )

    def fill(self, item: T) -> None:
        """
        Fill the grid with a given object.

        :param item:
            The item to fill the grid with.
        """
        self._data = [self.x_size * [item] for _ in range(self.y_size)]

    def map(self, func: Callable[[T], T]) -> None:
        for coord in self.coords:
            self[coord] = func(self[coord])

    def copy(self) -> "Grid":
        ret = type(self)(self.x_size, self.y_size)
        for coord in self.coords:
            ret[coord] = self[coord]
        return ret

    def deepcopy(self) -> "Grid":
        ret = type(self)(self.x_size, self.y_size)
        for coord in self.coords:
            ret[coord] = copy.deepcopy(self[coord])
        return ret

    def get_nbrs(self, coord: Coord, *, include_origin=False) -> Collection[Coord]:
        """
        Get a list of the coordinates of neighbouring cells.

        :param coord:
            The coordinate to get neighbours for.
        :param include_origin:
            Whether to include the original coordinate, coord, in the list.

        :return:
            List of coordinates within the boundaries of the grid.
        """
        x, y = coord
        nbrs = []
        for i in range(max(0, x - 1), min(self.x_size, x + 2)):
            for j in range(max(0, y - 1), min(self.y_size, y + 2)):
                nbrs.append((i, j))
        if not include_origin:
            nbrs.remove(coord)
        return nbrs


class Board(Grid[CellContents]):
    """
    Representation of a minesweeper board.

    Cells accessed by coordinate indexing, e.g. board[x, y], where (0, 0) is
    the top-left corner.

    May only contain instances of CellContents.
    """

    def __init__(self, x_size: int, y_size: int):
        """
        :param x_size:
            The number of columns.
        :param y_size:
            The number of rows.
        """
        super().__init__(x_size, y_size, fill=CellContents.Unclicked())

    def __str__(self):
        repr_grid = Grid.from_grid(self)
        repr_grid.map(lambda x: x.to_string())
        return str(repr_grid).replace("'", "")

    @classmethod
    def from_str(cls, board: str) -> "Board":
        """
        Create an instance from a string representation.

        :param board:
            The board string representation.
        :raise ValueError:
            If there is an invalid string representation for cell contents.
        :return:
            The created board instance.
        """
        grid = super().from_2d_iter(L.split() for L in board.splitlines() if L.strip())
        grid.map(CellContents.from_string)
        return cls.from_grid(grid)

    def reset(self):
        """Reset the board to the initial state."""
        self.fill(CellContents.Unclicked())
