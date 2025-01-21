#!/usr/bin/env python3
# A script to invoke mkdocs with the correct environment.
# Additionally supports deploying via mike:
#   ./mkdocs deploy

import errno
import os
import shutil
import sys
from pathlib import Path
from subprocess import call

contrib_doc_dir = Path(__file__).parent
root_dir = contrib_doc_dir.parents[1]
zstd_handler_dir = contrib_doc_dir / "mkdocstrings-zstd" / "src"
config_path = os.path.join(contrib_doc_dir, "mkdocs.yml")


build_dir = root_dir / "build-docs"

# Set PYTHONPATH for the mkdocstrings handler.
env = os.environ.copy()
path = env.get("PYTHONPATH")
env["PYTHONPATH"] = (path + ":" if path else "") + str(zstd_handler_dir)

args = sys.argv[1:]
args += ["-f", config_path]
sys.exit(call(["mkdocs"] + args, env=env))
