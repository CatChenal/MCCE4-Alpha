"""
CLI entry point for mcce4_gui.

Usage:
    python -m mcce4_gui                          # localhost:8080
    python -m mcce4_gui --port 9090              # custom port
    python -m mcce4_gui --host 0.0.0.0           # bind all interfaces
    python -m mcce4_gui --no-browser             # don't open browser
    python -m mcce4_gui --mcce-home ~/MCCE4      # set MCCE_HOME
"""

import argparse
import os
import sys
import threading
import webbrowser


def main():
    parser = argparse.ArgumentParser(
        prog="mcce4_gui",
        description="MCCE4 Web GUI — browser-based interface for MCCE4 pipeline",
    )
    parser.add_argument("--port", "-p", type=int, default=8080)
    parser.add_argument("--host", type=str, default="127.0.0.1")
    parser.add_argument("--no-browser", action="store_true")
    parser.add_argument("--mcce-home", type=str, default=None)
    parser.add_argument("--workdir", "-w", type=str, default=None)

    args = parser.parse_args()

    if args.mcce_home:
        os.environ["MCCE_HOME"] = args.mcce_home

    if args.workdir:
        os.makedirs(args.workdir, exist_ok=True)
        os.chdir(args.workdir)

    # Check deps
    try:
        import uvicorn
    except ImportError:
        print("Missing dependency: pip install uvicorn[standard]")
        sys.exit(1)

    try:
        import fastapi
    except ImportError:
        print("Missing dependency: pip install fastapi")
        sys.exit(1)

    url = f"http://{args.host}:{args.port}"
    is_ssh = bool(os.environ.get("SSH_CONNECTION"))

    print()
    print("  ╭───────────────────────────────────────────╮")
    print("  │            MCCE4 Web GUI v1.0.0            │")
    print("  ├───────────────────────────────────────────┤")
    print(f"  │  URL:      {url:<31s}│")
    if is_ssh or args.host == "127.0.0.1":
        print(f"  │  Tunnel:   ssh -L {args.port}:localhost:{args.port} user@host   │")
    print(f"  │  Workdir:  {os.getcwd()[:31]:<31s}│")
    mcce = os.environ.get("MCCE_HOME", args.mcce_home or "auto-detect")
    print(f"  │  MCCE:     {str(mcce)[:31]:<31s}│")
    print("  ╰───────────────────────────────────────────╯")
    print()

    if not args.no_browser and not is_ssh:
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()

    uvicorn.run(
        "mcce4_gui.app:app",
        host=args.host,
        port=args.port,
        log_level="info",
        reload=False,
    )


if __name__ == "__main__":
    main()
