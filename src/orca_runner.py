from pathlib import Path
import shutil
import subprocess


def find_orca() -> str:
    """Return absolute path to ORCA executable."""
    orca = shutil.which("orca")
    if orca is None:
        raise RuntimeError("ORCA executable not found in PATH")
    return orca


def run_orca(
    inp_file: Path,
    capture_out: bool = True,
    stream_output=None,
    stdout_file: Path | None = None
) -> dict:
    """
    Run ORCA (or GOAT) on the given input file.

    Args:
        inp_file: ORCA/GOAT input file
        capture_out: write a .out file (default True)
        stream_output: optional callable(line) for real-time output monitoring
        stdout_file: optional path to write stdout (default is inp_file.with_suffix('.out'))

    Returns:
        dict with keys: gbw, xyz, out
    """

    ORCA_CMD = find_orca()

    out_file = stdout_file or (inp_file.with_suffix(".out") if capture_out else None)

    if stream_output is None:
        # Buffered mode
        if out_file:
            with open(out_file, "w") as f:
                result = subprocess.run(
                    [ORCA_CMD, str(inp_file)],
                    stdout=f,
                    stderr=subprocess.STDOUT
                )
        else:
            result = subprocess.run([ORCA_CMD, str(inp_file)])

        if result.returncode != 0:
            raise RuntimeError(f"ORCA failed: check {out_file if out_file else 'terminal output'}")

    else:
        # Streaming mode
        if out_file is None:
            raise ValueError("stdout_file must be provided when using stream_output")

        with open(out_file, "w") as f:
            proc = subprocess.Popen(
                [ORCA_CMD, str(inp_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            for line in proc.stdout:
                f.write(line)
                stream_output(line.rstrip())

            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"ORCA failed: check {out_file}")

    return {
        "gbw": inp_file.with_suffix(".gbw"),
        "xyz": inp_file.with_suffix(".xyz"),
        "out": out_file
    }