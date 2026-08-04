from pathlib import Path
import os
import re
import signal
import shutil
import subprocess
import time
import queue
import threading
from collections.abc import Callable


def find_orca() -> str:
    """Return the ORCA executable from STORCA_ORCA_BIN or PATH."""
    configured = os.environ.get("STORCA_ORCA_BIN")
    if configured:
        executable = Path(configured).expanduser()
        if executable.is_file() and os.access(executable, os.X_OK):
            return str(executable.resolve())
        raise RuntimeError(f"STORCA_ORCA_BIN is not an executable file: {executable}")
    orca = shutil.which("orca")
    if orca is None:
        raise RuntimeError("ORCA executable not found in PATH")
    return orca


def validate_orca_output(out_file: Path) -> dict:
    """Require ORCA's own normal-termination marker, not just process status."""
    out_file = Path(out_file)
    if not out_file.is_file():
        raise RuntimeError(f"ORCA did not write an output file: {out_file}")
    text = out_file.read_text(errors="replace")
    normal = "ORCA TERMINATED NORMALLY" in text
    return {"normal_termination": normal, "output": str(out_file)}


def run_orca(
    inp_file: Path,
    capture_out: bool = True,
    stream_output: Callable[[str], None] | None = None,
    stdout_file: Path | None = None,
    extra_env: dict[str, str] | None = None,
    timeout_seconds: float | None = None,
    live_output: bool | None = None,
) -> dict:
    """
    Run ORCA (or GOAT) on the given input file.

    Args:
        inp_file: ORCA/GOAT input file
        capture_out: write a .out file (default True)
        stream_output: optional callable(line) for real-time output monitoring.
            A supplied callback takes precedence over ``live_output``.
        stdout_file: optional path to write stdout (default is inp_file.with_suffix('.out'))
        extra_env: environment variables required by this specific calculation
        live_output: mirror ORCA's native ``.out`` stream to this process's
            standard output.  Defaults to true; set ``STORCA_ORCA_LIVE_OUTPUT=0``
            to suppress it globally, or pass ``False`` for one job.

    Returns:
        dict with keys: gbw, xyz, out
    """

    inp_file = Path(inp_file).resolve()
    ORCA_CMD = find_orca()
    runtime_env = os.environ.copy()
    runtime_env.update(extra_env or {})

    out_file = stdout_file or (inp_file.with_suffix(".out") if capture_out else None)
    if live_output is None:
        live_output = os.environ.get("STORCA_ORCA_LIVE_OUTPUT", "1").strip().lower() not in {
            "0", "false", "no", "off",
        }

    # ORCA output is the most useful progress indicator for expensive jobs.
    # Retain the full file *and* mirror the exact stream to the calling
    # terminal by default.  The reader thread lets us continue enforcing a
    # timeout even when ORCA is temporarily quiet.
    if stream_output is None and live_output:
        stream_output = lambda line: print(line, end="", flush=True)

    if stream_output is None:
        # Buffered mode
        if out_file:
            with open(out_file, "w") as f:
                result = subprocess.run(
                    [ORCA_CMD, inp_file.name],
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    cwd=inp_file.parent,
                    env=runtime_env,
                    timeout=timeout_seconds,
                )
        else:
            result = subprocess.run(
                [ORCA_CMD, inp_file.name],
                cwd=inp_file.parent,
                env=runtime_env,
                timeout=timeout_seconds,
            )

        if result.returncode != 0:
            raise RuntimeError(f"ORCA failed: check {out_file if out_file else 'terminal output'}")
        if out_file and not validate_orca_output(out_file)["normal_termination"]:
            raise RuntimeError(f"ORCA exited without normal-termination marker: check {out_file}")

    else:
        # Streaming mode
        if out_file is None:
            raise ValueError("stdout_file must be provided when using stream_output")

        with open(out_file, "w") as f:
            proc = subprocess.Popen(
                [ORCA_CMD, inp_file.name],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=inp_file.parent,
                env=runtime_env,
            )

            assert proc.stdout is not None
            lines: queue.Queue[str | None] = queue.Queue()

            def _read_output() -> None:
                try:
                    for line in proc.stdout:
                        lines.put(line)
                finally:
                    lines.put(None)

            reader = threading.Thread(target=_read_output, name="storca-orca-output", daemon=True)
            reader.start()
            deadline = time.monotonic() + timeout_seconds if timeout_seconds is not None else None
            ended = False
            while not ended:
                try:
                    line = lines.get(timeout=0.2)
                except queue.Empty:
                    if deadline is not None and time.monotonic() > deadline:
                        proc.kill()
                        proc.wait()
                        raise subprocess.TimeoutExpired([ORCA_CMD, inp_file.name], timeout_seconds)
                    continue
                if line is None:
                    ended = True
                    continue
                f.write(line)
                f.flush()
                stream_output(line)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError(f"ORCA failed: check {out_file}")
        if not validate_orca_output(out_file)["normal_termination"]:
            raise RuntimeError(f"ORCA exited without normal-termination marker: check {out_file}")

    return {
        "gbw": inp_file.with_suffix(".gbw"),
        "xyz": inp_file.with_suffix(".xyz"),
        "out": out_file
    }
def run_rmg(
    input_file: Path,
    rmg_env: str | None = "rmg_env",
    rmg_command: str = "rmg.py",
    walltime: str | None = None,
    max_processes: int | None = None,
    max_iterations: int | None = None,
) -> dict:
    """
    Run RMG using 'conda run' so this works from any Python environment.

    Parameters
    ----------
    input_file : Path
        The RMG input script (e.g., 'input.py').
    rmg_env : str | None
        Name of the conda environment containing RMG. Set to ``None`` to run
        the configured RMG executable directly from PATH.
    rmg_command : str
        The RMG executable name (usually 'rmg.py').
    
    Raises
    ------
    subprocess.CalledProcessError
        If the RMG run fails.
    """
    # Ensure the working directory exists
    workdir = input_file.parent
    workdir.mkdir(parents=True, exist_ok=True)

    if walltime is not None and not re.fullmatch(r"\d{2,}:\d{2}:\d{2}:\d{2}", walltime):
        raise ValueError("RMG walltime must use DD:HH:MM:SS, for example 00:00:10:00")
    if max_processes is not None and max_processes < 1:
        raise ValueError("RMG max_processes must be at least one")
    if max_iterations is not None and max_iterations < 1:
        raise ValueError("RMG max_iterations must be at least one")
    options: list[str] = []
    if walltime:
        options.extend(["-t", walltime])
    if max_processes:
        options.extend(["-n", str(max_processes)])
    if max_iterations:
        options.extend(["-i", str(max_iterations)])
    if rmg_env:
        cmd = ["conda", "run", "-n", rmg_env, rmg_command, *options, str(input_file.name)]
    else:
        cmd = [rmg_command, *options, str(input_file.name)]

    # Run RMG in the directory containing input_file
    result = subprocess.run(
        cmd,
        cwd=workdir,           # <-- ensures RMG.log is created here
        capture_output=True,   # capture stdout/stderr
        text=True,
    )
    if result.returncode != 0:
        error = RuntimeError(f"RMG failed with exit code {result.returncode}")
        error.stdout = result.stdout  # type: ignore[attr-defined]
        error.stderr = result.stderr  # type: ignore[attr-defined]
        raise error
    return {"command": cmd, "stdout": result.stdout, "stderr": result.stderr}


def rmg_walltime_seconds(walltime: str | None) -> float | None:
    """Convert RMG's DD:HH:MM:SS walltime form to seconds."""
    if walltime is None:
        return None
    days, hours, minutes, seconds = (int(part) for part in walltime.split(":"))
    return float((((days * 24) + hours) * 60 + minutes) * 60 + seconds)


def run_rmg_supervised(
    input_file: Path,
    rmg_env: str | None = "rmg_env",
    rmg_command: str = "rmg.py",
    walltime: str | None = None,
    max_processes: int | None = None,
    max_iterations: int | None = None,
    *,
    hard_timeout_seconds: float | None = None,
    heartbeat_timeout_seconds: float | None = None,
    poll_seconds: float = 2.0,
    start_method: str = "auto",
    progress_interval_seconds: float = 30.0,
) -> dict:
    """Run RMG under an external wall-clock/heartbeat supervisor.

    RMG's own ``-t`` limit controls model generation.  This wrapper protects
    the calling workflow from a wedged executable or conda wrapper.  It does
    not infer convergence; callers must inspect retained RMG artifacts.
    """
    input_file = Path(input_file)
    if walltime is not None and not re.fullmatch(r"\d{2,}:\d{2}:\d{2}:\d{2}", walltime):
        raise ValueError("RMG walltime must use DD:HH:MM:SS")
    options: list[str] = []
    if walltime:
        options.extend(["-t", walltime])
    if max_processes:
        options.extend(["-n", str(max_processes)])
    if max_iterations:
        options.extend(["-i", str(max_iterations)])
    if start_method not in {"auto", "fork", "spawn"}:
        raise ValueError("RMG start_method must be auto, fork, or spawn")
    launcher = Path(__file__).resolve().parent.parent / "storca" / "rmg_launcher.py"
    if rmg_env and launcher.is_file():
        command = ["conda", "run", "-n", rmg_env, "python", str(launcher),
                   "--storca-start-method", start_method, *options, input_file.name]
    else:
        command = [rmg_command, *options, input_file.name]
    limit = hard_timeout_seconds
    if limit is None:
        native = rmg_walltime_seconds(walltime)
        limit = native + 120.0 if native is not None else None
    stdout_path, stderr_path = input_file.parent / "rmg.supervisor.stdout.log", input_file.parent / "rmg.supervisor.stderr.log"
    started = time.monotonic()
    heartbeat = input_file.parent / "RMG.log"
    last_change = started
    last_mtime = heartbeat.stat().st_mtime if heartbeat.is_file() else None
    stop_reason = None
    next_progress = started
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(command, cwd=input_file.parent, stdout=stdout, stderr=stderr,
                                   start_new_session=True, text=True)
        while process.poll() is None:
            now = time.monotonic()
            mtime = heartbeat.stat().st_mtime if heartbeat.is_file() else None
            if mtime != last_mtime:
                last_mtime, last_change = mtime, now
            if limit is not None and now - started > limit:
                stop_reason = "hard_timeout"
            elif heartbeat_timeout_seconds is not None and now - last_change > heartbeat_timeout_seconds:
                stop_reason = "heartbeat_timeout"
            if stop_reason:
                os.killpg(process.pid, signal.SIGINT)
                try:
                    process.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                break
            if now >= next_progress:
                elapsed = int(now - started)
                state = "waiting for RMG.log" if not heartbeat.is_file() else "RMG.log active"
                print(f"[STORCA {elapsed // 60:02d}:{elapsed % 60:02d}] RMG model generation running ({state})", flush=True)
                next_progress = now + max(5.0, progress_interval_seconds)
            time.sleep(max(0.1, poll_seconds))
    return {
        "command": command, "stdout": stdout_path.read_text(errors="replace"),
        "stderr": stderr_path.read_text(errors="replace"), "returncode": process.returncode,
        "supervisor": {"status": "interrupted" if stop_reason else "completed",
                       "stop_reason": stop_reason, "elapsed_seconds": time.monotonic() - started,
                       "hard_timeout_seconds": limit, "heartbeat_timeout_seconds": heartbeat_timeout_seconds},
    }


def run_arkane(
    input_file: Path,
    rmg_env: str | None = "rmg_env",
    arkane_command: str = "Arkane.py",
) -> dict:
    """Run Arkane in the RMG environment and retain its console streams."""
    input_file = Path(input_file)
    input_file.parent.mkdir(parents=True, exist_ok=True)
    command = (["conda", "run", "-n", rmg_env, arkane_command, input_file.name]
               if rmg_env else [arkane_command, input_file.name])
    result = subprocess.run(command, cwd=input_file.parent, capture_output=True, text=True)
    if result.returncode != 0:
        error = RuntimeError(f"Arkane failed with exit code {result.returncode}")
        error.stdout = result.stdout  # type: ignore[attr-defined]
        error.stderr = result.stderr  # type: ignore[attr-defined]
        raise error
    return {"command": command, "stdout": result.stdout, "stderr": result.stderr}
