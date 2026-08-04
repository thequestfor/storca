#!/usr/bin/env python
"""Launch RMG with a macOS-safe multiprocessing policy.

RMG initializes its database in the parent process.  Python's macOS ``spawn``
default starts clean workers which therefore lack that database; Unix ``fork``
preserves the initialized state and is RMG's normal Linux behaviour.
"""

from __future__ import annotations

import multiprocessing as multiprocessing
import sys


def main() -> None:
    requested = "auto"
    if "--storca-start-method" in sys.argv:
        index = sys.argv.index("--storca-start-method")
        try:
            requested = sys.argv[index + 1]
        except IndexError as error:
            raise SystemExit("--storca-start-method requires a value") from error
        del sys.argv[index:index + 2]
    if requested not in {"auto", "fork", "spawn"}:
        raise SystemExit("STORCA RMG start method must be auto, fork, or spawn")
    method = "fork" if requested == "auto" and sys.platform == "darwin" else requested
    if method != "auto":
        multiprocessing.set_start_method(method)
    # Import only after selecting the process start method.
    from rmgpy.__main__ import main as rmg_main
    rmg_main()


if __name__ == "__main__":
    main()
