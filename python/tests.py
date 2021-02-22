import pytest

from .utils import Board, CellContents, Grid


@pytest.mark.parametrize(
    "test_input, expected",
    [
        ("#", CellContents.Unclicked()),
        (".", CellContents.Num(0)),
        ("0", CellContents.Num(0)),
        ("1", CellContents.Num(1)),
        ("5", CellContents.Num(5)),
        ("F", CellContents.Flag(1)),
        ("F1", CellContents.Flag(1)),
        ("F2", CellContents.Flag(2)),
        ("X2", CellContents.WrongFlag(2)),
        ("M1", CellContents.Mine(1)),
        ("!", CellContents.HitMine(1)),
    ],
)
def test_cell_contents_from_str(test_input: str, expected: CellContents):
    assert CellContents.from_string(test_input) == expected


@pytest.mark.parametrize(
    "test_input, expected",
    [
        (CellContents.Unclicked(), "#"),
        (CellContents.Num(0), "."),
        (CellContents.Num(1), "1"),
        (CellContents.Num(5), "5"),
        (CellContents.Flag(1), "F1"),
        (CellContents.Flag(2), "F2"),
        (CellContents.WrongFlag(2), "X2"),
        (CellContents.Mine(1), "M1"),
        (CellContents.HitMine(1), "!1"),
    ],
)
def test_cell_contents_to_str(test_input: CellContents, expected: str):
    assert test_input.to_string() == expected


def test_create_grid():
    Grid(3, 4)


def test_create_board():
    Board(3, 4)


def test_create_board_from_str():
    board = Board.from_str(
        """
            #  . 1   2
            
            #  3 F   !2
        """
    )
    assert board.x_size == 4
    assert board.y_size == 2
    assert list(board.values) == [
        CellContents.Unclicked(),
        CellContents.Num(0),
        CellContents.Num(1),
        CellContents.Num(2),
        CellContents.Unclicked(),
        CellContents.Num(3),
        CellContents.Flag(1),
        CellContents.HitMine(2),
    ]
