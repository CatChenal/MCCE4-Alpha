# mcce4_gui — Web GUI for MCCE4

A modern, browser-based interface for the MCCE4 protein electrostatics pipeline. Runs as a lightweight web server on the compute node and streams everything to your browser — no X11 forwarding needed.

## Installation

The GUI lives inside the MCCE4-Alpha repository as a subpackage:

```
MCCE4-Alpha/
├── bin/
│   ├── mcce4-gui          ← launcher script
│   ├── mcce4_gui/         ← this package
│   │   ├── __init__.py
│   │   ├── __main__.py
│   │   ├── app.py
│   │   ├── config.py
│   │   ├── pipeline.py
│   │   ├── slurm.py
│   │   ├── analysis.py
│   │   └── static/
│   │       └── index.html
│   ├── step1.py
│   ├── step2.py
│   └── ...
```

Install the Python dependencies:

```bash
pip install fastapi uvicorn[standard] python-multipart websockets
```

## Quick Start

```bash
# From any MCCE4 working directory:
cd /path/to/protein_workdir
mcce4-gui --no-browser --mcce-home ~/MCCE4
```

Then SSH tunnel from your local machine and open `localhost:8080`:

```bash
ssh -L 8080:localhost:8080 user@server
# open http://localhost:8080
```

## Features

### Pipeline Tab
- **Visual step tracker** — Preparation → Rotamers → Energies → Monte Carlo with live progress
- **3D structure preview** — Interactive protein visualization via 3Dmol.js (no PyMOL needed)
- **Server file browser** — Browse and select PDB files without typing paths
- **Live log streaming** — Real-time output via WebSocket
- **Stop & resume** — Interrupt long computations and pick up where you left off

### Configuration Tab
- **SLURM submission** — Generate sbatch scripts, submit jobs, monitor status
- **Script preview** — Review the exact sbatch script before submitting
- **Persistent configs** — Save/load parameter sets as JSON

### Analysis Tab
- **Titration curves** — Interactive Plotly charts of occupancy vs pH
- **Net charge plot** — Total charge vs pH from sum_crg.out
- **pKa table** — Sortable results with residue names, pKa values, Hill coefficients
- **Output file checker** — Verify which MCCE4 output files exist

### Light & Dark Mode
Automatically follows your OS preference. The 3D viewer, Plotly charts, and all UI elements adapt.

## CLI Options

```
mcce4-gui [OPTIONS]

  --port, -p INT       Port (default: 8080)
  --host TEXT          Host (default: 127.0.0.1)
  --no-browser         Don't auto-open browser
  --mcce-home TEXT     MCCE4 installation path
  --workdir, -w TEXT   Working directory
```

## REST API

All functionality is accessible programmatically:

| Endpoint | Method | Description |
|---|---|---|
| `/api/config` | GET/POST | Read/update configuration |
| `/api/config/save` | POST | Persist to disk |
| `/api/pipeline/run` | POST | Start pipeline |
| `/api/pipeline/stop` | POST | Stop pipeline |
| `/api/pipeline/resume` | POST | Resume pipeline |
| `/api/pipeline/status` | GET | Get step statuses |
| `/api/files/browse` | GET | Browse server filesystem |
| `/api/upload/pdb` | POST | Upload PDB file |
| `/api/pdb/content` | GET | Get PDB content for viewer |
| `/api/slurm/submit` | POST | Submit SLURM job |
| `/api/slurm/status` | GET | Poll job status |
| `/api/slurm/preview` | GET | Preview sbatch script |
| `/api/slurm/cancel` | POST | Cancel SLURM job |
| `/api/analysis/pka` | GET | Parsed pK.out data |
| `/api/analysis/sum_crg` | GET | Parsed charge data |
| `/api/analysis/conformers` | GET | Parsed head3.lst |
| `/api/analysis/available` | GET | Check output files |
| `/ws` | WebSocket | Real-time log stream |

## Requirements

- Python 3.8+
- FastAPI, Uvicorn, python-multipart, websockets
- MCCE4 installed and accessible
- A modern web browser

## License

MIT — Same as MCCE4-Alpha
