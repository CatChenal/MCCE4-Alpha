#!/usr/bin/env python
"""
Module: downloads.py

Provides functions to download files with the requests library.
"""
import argparse
import gzip
import logging
from pathlib import Path
import requests
import shutil
from typing import Tuple, Union


logging.basicConfig(format="[ %(levelname)s ] - %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# reset libs logging to higher level so that no unnecessary message is logged
# when this module is used.
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
SMILES_QRY = """
query molecules($ids: [String!]!) {
  chem_comps(comp_ids: $ids) {
    chem_comp {
      id
      name
      formula
      pdbx_formal_charge
    }
    rcsb_chem_comp_descriptor {
      SMILES
    }
  }
}
"""


def get_rcsb_lig_smiles(ligands:list, output_dir:str=Path("./prerun")):
    """Get the smiles strings for ligands (3-letter code).
    Save the list of [id, ligand] to csv file in the output_dir
    as smiles.csv.
    """
    variables = {"ids": ligands}
    payload = {
        "query": SMILES_QRY,
        "variables": variables
    }
    response = requests.post(RCSB_GRAPHQL, json=payload)
    if not response.ok:
        print(f"Failed smiles request: {response.reason}")
        return
    data = response.json()
    # concate long string
    out = [(f'{entry["chem_comp"]["id"]},'
            f'{entry["rcsb_chem_comp_descriptor"]["SMILES"]},'
            f'{entry["chem_comp"]["formula"]}\n')
           for entry in data["data"]["chem_comps"]]
    if out:
        with open(output_dir.joinpath("smiles.csv"), "w") as fo:
            fo.writelines(out)
    else:
        print("No smiles data")
    return


def rcsb_download(pdb_fname: str) -> requests.Response:
    url_rscb = "https://files.rcsb.org/download/" + pdb_fname
    return requests.get(url_rscb, allow_redirects=True)


def rcsb_download_fasta(pdb_id: str) -> requests.Response:
    fasta_url = "https://www.rcsb.org/fasta/entry/" + pdb_id 
    return requests.get(fasta_url, allow_redirects=True)


def get_rcsb_fasta(pdbid: str, output_dir: str = Path.cwd(), verbose: bool=False) -> int:
    """Given a pdb id, download its fasta file if not found."""
    pdbid = pdbid.upper()
    pdb = Path(output_dir).joinpath(f"{pdbid}.fasta")
    if pdb.exists():
        if verbose:
            print(f" {pdbid} fasta file already exists.")
        return 0

    response = rcsb_download_fasta(pdbid)
    if response.status_code < 400:
        with open(pdb, "w") as fo:
            fo.write(response.text)
        if verbose:
            print(f" Saved {pdbid} fasta file.")
        return 0
    else:
        print(f" ERROR: Could not download fasta file of {pdbid}:\n{response.reason}.")
        return 1

    
def get_rcsb_pdb(pdbid: str,
                 get_bioassembly: bool = True,
                 bioassembly_id: int = 1,
                 cif_format: bool = False) -> Union[Path, Tuple[None, str]]:
    """Given a pdb id, possibly download the pdb file containing the biological
    assembly with given bioassembly_id from rcsb.org.
    If cif_format is True, the download of the bioassembly in .cif format is attempted,
    otherwise, the file download is attempted in this order:
     If get_bioassembly is False or the download of the bioassembly in pdb format fails,
     then the download of the standard pdb file is attempted.

    Arguments:
     - pdbid (str): The pdb id to download
     - get_bioassembly (bool, True): Whether to attempt the bioassembly or the full structure
       download
     - bioassembly_id (int, 1): Which bioassembly to download if get_bioassembly is True.
     - cif_format (bool, False): Whether to download the file in .cif format directly.

    Returns:
     - The path to the downloaded pdb file or a tuple signifying and error occured with values
       (None, error_message)
    """
    pdbid = pdbid.lower()

    if cif_format:
        # only get the cif format with or w/o bioassembly
        if not get_bioassembly:
            pdb = f"{pdbid}.cif.gz"
        else:
            pdb = f"{pdbid}-assembly{bioassembly_id}.cif.gz"

        r0 = rcsb_download(pdb)
        if r0.status_code < 400:
            with open(pdb, "wb") as fo:
               fo.write(r0.content)
            decomp = decompress_gz(Path(pdb))
            Path(pdb).unlink()
            return Path(decomp).resolve()
        else:
            logger.warning(f"Could not download the structure in cif format: {r0.reason}.")
            return None, "Error: Could not download the structure in cif format: {r0.reason}."
    else:
        # try pdb, then .cif format as a fallback
        if not get_bioassembly:
            pdb_names = pdbid + ".pdb", f"{pdbid}.cif.gz"
        else:
            pdb_names = pdbid + ".pdb1", f"{pdbid}-assembly{bioassembly_id}.cif.gz"
        which_kinds = [False, False]
    
        # first attempt: pdb format
        pdb_file = pdb_names[0]
        r0 = rcsb_download(pdb_file)
        if r0.status_code < 400:
            which_kinds[0] = True
            with open(pdb_file, "wb") as fo:
                fo.write(r0.content)
        else:
            logger.warning(f"Could not download the standard pdb: {r0.reason}. Trying cif.")
            
        if not which_kinds[0]:
            # second attempt: cif format
            pdb_file = pdb_names[1]
            r1 = rcsb_download(pdb_file)
            if r1.status_code < 400:
                which_kinds[1] = True
                with open(pdb_file, "wb") as fo:
                    fo.write(r1.content)
            else:
                logger.warning(f"Could not download the structure in cif format: {r1.reason}.")
                
            if not which_kinds[1]:  # both False
                if get_bioassembly:
                    msg = """"Could neither download the bioassembly structure in pdb nor in cif format. 
                    Perhaps, try again with: -get_bioassembly False."""
                    logger.error(msg)
                else:  
                    logger.error("Could neither download the structure in pdb nor in cif format.")
                return None, "Error: Could neither download the structure in pdb nor in cif format."
            else:
                decomp = decompress_gz(Path(pdb_file))
                Path(pdb_file).unlink()
                return Path(decomp).resolve()
        else:
            if not get_bioassembly:
                return Path(pdb_file).resolve()
            else:
                pdb = pdbid + ".pdb"
                shutil.move(pdb_file, pdb)
                return Path(pdb).resolve()


def opt2bool(opt) -> bool:
    if opt =="True":
        return True
    return False
    
def cli_parser():
    p = argparse.ArgumentParser(prog="getpdb",
                                description=("Download one or more pdb files (bioassembly by default)"
                                             " from the RSCB download service."))
    p.add_argument(
        "pdbid",
        nargs="+",
        default=[],
        help="Specify the pdb ID(s), e.g.: 1ots 4lzt 1FAT",
    )
    # type cannot be bool to be updatable
    p.add_argument(
        "-get_bioassembly",
        default="True",
        type=opt2bool,
        help="Whether to attempt the bioassembly (default) or the full pdb download"
    )
    p.add_argument(
        "-bioassembly_id",
        type=int,
        default=1,
        help="Which bioassembly to download (default: 1)."
    )
    p.add_argument(
        "--cif",
        default=False,
        action="store_true",
        help="Whether to downloaded a bioassembly file in cif format (default: False)"
    )

    return p


def getpdb_cli(argv=None):
    """Cli function for the `getpdb` tool.
    """
    p = cli_parser()
    args = p.parse_args(argv)

    for pdbid in [id.lower() for id in args.pdbid]:
        out = get_rcsb_pdb(pdbid,
                           get_bioassembly=args.get_bioassembly,
                           bioassembly_id=args.bioassembly_id,
                           cif_format=args.cif)
        if isinstance(out, tuple):
            logger.info(f"Download failed: {pdbid}; {out[1]}")
        else:
            logger.info(f"Download completed: {out.name}")


def fasta_parser():
    p = argparse.ArgumentParser(prog="getfasta",
                                description=("Download one or more fasta files from RSCB."))
    p.add_argument(
        "pdbid",
        nargs="+",
        default=[],
        help="Specify the pdb ID(s), e.g.: 1ots 4lzt 1FAT",
    )
    return p


def getfasta_cli(argv=None):
    """Cli function for the `getfasta` tool.
    """
    p = fasta_parser()
    args = p.parse_args(argv)

    for pdbid in [id.upper() for id in args.pdbid]:
        out = get_rcsb_fasta(pdbid)
        if out:
            logger.info(f"Download failed: {pdbid}")


def decompress_gz(gfp: Path) -> Path:
    """Decompressed the gzipped file given by its filepath, gfp.
    """
    fpo = gfp.parent.joinpath(f"{gfp.stem}")
    with gzip.open(gfp, "rb") as f_in:
        with open(fpo, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    return fpo
