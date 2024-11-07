# Copyright 2023 Michael Davidson (UCSD), Xiqiang Liu (UCSD)

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
logger.setLevel(_logging.INFO)
logger.propagate = False
if not logger.hasHandlers():
    ch = _logging.StreamHandler()
    ch.setFormatter(
        CustomFormatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )
    logger.addHandler(ch)
