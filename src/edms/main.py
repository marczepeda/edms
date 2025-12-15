#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
'''
Module: main.py
Author: Marc Zepeda
Created: 2025-04-12
Description: Endogenous Deep Mutational Scans

Usage:
[Show post-install notice only once per version]
- try: importlib.metadata to get current version
- _user_config_dir(appname): Get cross-platform user config directory
- _maybe_show_first_run_notice(appname): Show a post-install notice only once per version

[Main method]
- main(): Endogenous Deep Mutational Scans
'''
# Import packages
from __future__ import annotations # NEEDS TO BE FIRST LINE
import json, os, sys
from pathlib import Path
import argparse
import argcomplete
import ast
import datetime
from rich_argparse import RichHelpFormatter
from rich import print as rprint
import sys

from . import config
from . import utils
from .bio import cli as bio_cli
from .data import cli as data_cli
from .gen import cli as gen_cli

# Show post-install notice only once per version
try:
    # Py3.8+: importlib.metadata in stdlib
    from importlib.metadata import version, PackageNotFoundError
except Exception:  # pragma: no cover
    version = lambda _: "0.0.0"  # fallback

def _user_config_dir(appname: str) -> Path:
    """
    _user_config_dir(): Cross-platform config dir without extra deps.
    
    Parameters:
    appname (str): Name of the application (used to create app-specific subdirectory)
    """
    if os.name == "nt":  # Windows
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
    else:  # POSIX
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    return base / appname

def _maybe_show_first_run_notice(appname: str = "edms") -> None:
    """
    _maybe_show_first_run_notice(): Show a post-install notice only once per version.
    
    Parameters:
    appname (str): Name of the application (used to create app-specific subdirectory)
    """
    # Allow users/CI to suppress the message
    if os.environ.get("EDMS_NO_FIRST_RUN") == "1":
        return

    cfg_dir = _user_config_dir(appname)
    state_file = cfg_dir / "state.json"
    cfg_dir.mkdir(parents=True, exist_ok=True)

    # Determine current installed version (best effort)
    try:
        cur_ver = version("edms")
    except Exception:
        cur_ver = None

    state = {}
    if state_file.exists():
        try:
            state = json.loads(state_file.read_text() or "{}")
        except json.JSONDecodeError:
            state = {}

    last_shown_for = state.get("first_run_notice_shown_for")

    # Show only if never shown for this version (or ever)
    if cur_ver and last_shown_for == cur_ver:
        return

    # ---- Your message goes here ----
    rprint("[green]\nRecommended post-install:\n\n  1) edms autocomplete\n     Follow CLI instructions\n\n  2) conda create -n umi_tools umi_tools cutadapt samtools bowtie2 fgbio\n     Required for edms fastq {extract_umis,trim_motifs,make_sams,make_bams,bam_umi_tags,group_umis,consensus_umis,bam_to_fastq}\n[/green]",
           file=sys.stderr,
    )
    # --------------------------------

    state["first_run_notice_shown_for"] = cur_ver or "unknown"
    try:
        state_file.write_text(json.dumps(state, indent=2))
    except OSError:
        # Non-fatal: if we can't write, just skip persisting
        pass

# Main method
def main():
    '''
    main(): Endogenous Deep Mutational Scans
    '''
    rprint("[green]Project: Endogenous Deep Mutational Scans (EDMS)[/green]")
    _maybe_show_first_run_notice()

    # Custom formatter for rich help messages
    class MyFormatter(RichHelpFormatter):
        styles = {
            "argparse.prog": "green",           # program name
            "argparse.args": "cyan",            # positional arguments
            "argparse.option": "",              # options like --flag
            "argparse.metavar": "dark_magenta", # meta variable (actual function argument name)
            "argparse.help": "blue",            # help text
            "argparse.text": "green",           # normal text in help message
            "argparse.groups": "red",           # group titles
            "argparse.description": "",         # description at the top
            "argparse.epilog": "",              # ... -h; epilog at the bottom
            "argparse.syntax": "white",         # []
        }

    # Add parser and subparsers
    parser = argparse.ArgumentParser(prog="edms", description="Endogenous Deep Mutational Scans (EDMS)", formatter_class=MyFormatter)
    subparsers = parser.add_subparsers(dest="command") # dest="command" required for autocomplete

    '''
    edms.utils:
    - run_bundled_script(): Run a bundled script from the package script resources.
    '''
    subparsers.add_parser("autocomplete", help="Enable edms autocomplete", description="Enable edms autocomplete", formatter_class=MyFormatter) # Run autocomplete script (args.command == "autocomplete")

    '''
    edms.config:
    - get_info: Retrieve information based on id
    - set_info: Set information based on id
    - del_info: Delete information based on id
    '''
    parser_config = subparsers.add_parser("config", help="Configuration", description="Configuration", formatter_class=MyFormatter)
    subparsers_config = parser_config.add_subparsers()

    # Add subparsers for config module
    parser_config_get_info = subparsers_config.add_parser("get", help="Retrieve information based on id", description="Retrieve information based on id", formatter_class=MyFormatter)
    parser_config_set_info = subparsers_config.add_parser("set", help="Set information based on id", description="Set information based on id", formatter_class=MyFormatter)
    parser_config_del_info = subparsers_config.add_parser("del", help="Delete information based on id", description="Delete information based on id", formatter_class=MyFormatter)

    # get_info() arguments
    parser_config_get_info.add_argument("--id", type=str, help="Identifier from/for configuration file (e.g., SEER_API_KEY or USPTO_API_KEY)")
    
    # Commond set_info() and del_info() arguments
    for parser_common in [parser_config_set_info, parser_config_del_info]:
        parser_common.add_argument("--id", type=str, help="Identifier from/for configuration file (e.g., SEER_API_KEY or USPTO_API_KEY)", required=True)
    
    # set_info() arguments
    parser_config_set_info.add_argument("--info", type=str, help="Information for configuration file (str, dict, list, set, or tuple; e.g., your_api_key)", required=True)
    
    # default functions
    parser_config_get_info.set_defaults(func=config.get_info)
    parser_config_set_info.set_defaults(func=config.set_info)
    parser_config_del_info.set_defaults(func=config.del_info)

    # General module subparser
    gen_cli.add_subparser(subparsers, MyFormatter)

    # Data module subparser
    data_cli.add_subparser(subparsers, MyFormatter)

    # Bio module subparser
    bio_cli.add_subparser(subparsers, MyFormatter)

    # Enable autocomplete
    argcomplete.autocomplete(parser)

    # Parse all arguments
    args = parser.parse_args()
    if args.command == 'autocomplete': # Run autocomplete script
        sys.exit(utils.run_bundled_script())
    args_dict = vars(args)
    func = args_dict.pop("func")
    args_dict.pop("command", None)  # Remove 'command' if it exists (required for autocomplete)
    func(**args_dict)