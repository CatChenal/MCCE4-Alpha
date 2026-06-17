#!/usr/bin/env python3

"""
Module: launch_gui.py

"""
import argparse
import os
from pathlib import Path
from re import search as re_search
import signal
import socket
import subprocess
import sys
from typing import Union

try:
    import streamlit as st
except ImportError:
    from mcce4.mcce_ftpl_agent import gui_reqs_cmd
    print("\nMissing GUI dependencies!\n",gui_reqs_cmd())
    sys.exit(1)


GUI_APP_FP = Path(__file__).parent.joinpath("gui", "app.py")


def _is_port_in_use(port: int) -> bool:
    """Check if a TCP port is in use."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("localhost", port)) == 0


def _find_pid_on_port(port: int):
    """Find PID of process using the given port. Returns int or None."""
    try:
        result = subprocess.run(
            ["lsof", "-ti", f":{port}"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            # May return multiple PIDs; take the first
            return int(result.stdout.strip().split()[0])
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        pass

    # Fallback: try ss
    try:
        result = subprocess.run(
            ["ss", "-tlnp", f"sport = :{port}"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            m = re_search(r'pid=(\d+)', result.stdout)
            if m:
                return int(m.group(1))
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    return None


def launch_gui(args: Union[argparse.Namespace, dict]):
    """Launch the Streamlit web GUI.
    """
    if not GUI_APP_FP.exists():
        print(f"ERROR: GUI app not found at {GUI_APP_FP!s}")
        sys.exit(1)
    gui_path = str(GUI_APP_FP)

    # Check if port is already in use and handle it
    port = 8501
    if _is_port_in_use(port):
        print(f"⚠  Port {port} is already in use (previous Streamlit session?).")
        print()

        # Try to find the PID
        pid = _find_pid_on_port(port)
        if pid:
            print(f"   Found process PID {pid} on port {port}.")
            try:
                response = input("   Kill it and restart? [Y/n] ").strip().lower()
            except EOFError:
                response = "y"

            if response in ("", "y", "yes"):
                try:
                    os.kill(pid, signal.SIGTERM)
                    print(f"   Killed PID {pid}. Waiting for port to free...")
                    import time
                    for _ in range(10):
                        time.sleep(0.5)
                        if not _is_port_in_use(port):
                            break
                    if _is_port_in_use(port):
                        os.kill(pid, signal.SIGKILL)
                        time.sleep(1)
                except ProcessLookupError:
                    pass  # Already dead
                except PermissionError:
                    print("   Permission denied. Run manually:")
                    print(f"     kill {pid}")
                    print("   Then rerun this command.")
                    sys.exit(1)
            else:
                print("\n   To free the port manually, run:")
                print(f"     kill {pid}")
                print("   Or use a different port:")
                print(f"     streamlit run {gui_path} --server.port 8502 -- {os.path.abspath(args.pdb)}")
                sys.exit(1)
        else:
            print(f"   Could not identify the process. To free port {port}, run:")
            print("     pkill -f 'streamlit run'")
            print(f"     # or: kill $(lsof -ti:{port})")
            print("   Then rerun this command.")
            sys.exit(1)

    if isinstance(args, dict):
        args = argparse.Namespace(**args)

    # Pass PDB path via environment so Streamlit can access it
    os.environ["MCCE_AGENT_PDB"] = os.path.abspath(args.pdb)
    os.environ["MCCE_AGENT_PH"] = str(args.ph)
    os.environ["MCCE_AGENT_CHARGE_METHOD"] = args.charge_method

    print("🌐 Launching MCCE4 Topology Agent GUI...")
    print(f"   PDB: {args.pdb}")
    print(f"   Open your browser to: http://localhost:{port}")
    print(f"   (If remote: ssh -L {port}:localhost:{port} user@server)\n")

    try:
        subprocess.run([
            sys.executable, "-m", "streamlit", "run", gui_path,
            "--server.headless", "true",
            "--server.port", str(port),
            "--browser.gatherUsageStats", "false",
        ])
    except KeyboardInterrupt:
        print("\n  GUI stopped.")
    except FileNotFoundError:
        print("ERROR: Streamlit not installed — install: conda install -c conda-forge streamlit")
        sys.exit(1)
