# Copyright 2023, 2025 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import contextlib
import logging as _logging

color2num = {
    "gray": 30,
    "red": 31,
    "green": 32,
    "yellow": 33,
    "blue": 34,
    "magenta": 35,
    "cyan": 36,
    "white": 37,
    "crimson": 38,
}


def colorize(
    string: str, color: str, bold: bool = False, highlight: bool = False
) -> str:
    """Returns string surrounded by appropriate terminal colour codes to print colorized text.

    Args:
        string: The message to colourise
        color: Literal values are gray, red, green, yellow, blue, magenta, cyan, white, crimson
        bold: If to bold the string
        highlight: If to highlight the string

    Returns:
        Colourised string
    """
    attr = []
    num = color2num[color]
    if highlight:
        num += 10
    attr.append(str(num))
    if bold:
        attr.append("1")
    attrs = ";".join(attr)
    return f"\x1b[{attrs}m{string}\x1b[0m"


class StdoutToLoggerRedirect:
    """A class to redirect stdout to a logger.

    Args:
        logger: The logger to redirect stdout to. If not provided, a new logger
            will be created.
        level: The logging level to use when redirecting stdout.
            Default is logging.INFO.

    Example:
    >>> import logging
    >>> with contextlib.redirect_stdout(StdoutToLoggerRedirect()):
    ...     print("Hello, world!")
    """

    def __init__(
        self, logger: _logging.Logger | None = None, level: int = _logging.INFO
    ):
        self.logger = logger or _logging.getLogger(__name__)
        self.level = level

    def write(self, msg: str):
        if msg and not msg.isspace():
            self.logger.log(self.level, msg)

    def flush(self):
        pass


@contextlib.contextmanager
def redirect_stdout_to_logger(logger=None, level=_logging.INFO):
    """Context manager to redirect stdout to a logger.

    Args:
        logger: The logger to redirect stdout to. If not provided, a new logger
            will be created.
        level: The logging level to use. Default is logging.INFO.

    Example:
    >>> import logging
    >>> with redirect_stdout_to_logger():
    ...     print("Hello, world!")
    """
    with contextlib.redirect_stdout(StdoutToLoggerRedirect(logger, level)):
        yield


class CustomFormatter(_logging.Formatter):
    # https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output

    def format(self, record):
        s = super().format(record)
        if record.levelno == _logging.WARNING:
            s = colorize(s, "yellow", True)
        elif record.levelno == _logging.ERROR:
            s = colorize(s, "red", True, True)
        elif record.levelno == _logging.INFO:
            s = colorize(s, "green")
        elif record.levelno == _logging.DEBUG:
            s = colorize(s, "blue")
        return s


logger = _logging.getLogger("geodata")

# Only set the level if it hasn't been set yet
if logger.level == _logging.NOTSET:
    logger.setLevel(_logging.INFO)

logger.propagate = False
if not logger.hasHandlers():
    ch = _logging.StreamHandler()
    ch.setFormatter(
        CustomFormatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(ch)

__all__ = ["logger", "redirect_stdout_to_logger"]
