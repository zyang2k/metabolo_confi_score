"""Wrapper: reads token from /tmp/.mw_token and runs refetch."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

token = open("/tmp/.mw_token").read().strip()
dry = "--dry-run" in sys.argv
limit = None
for i, a in enumerate(sys.argv):
    if a == "--limit" and i + 1 < len(sys.argv):
        limit = int(sys.argv[i + 1])

# Patch argv so argparse sees the right values
new_argv = [sys.argv[0], "--token", token]
if dry:
    new_argv.append("--dry-run")
if limit is not None:
    new_argv += ["--limit", str(limit)]
sys.argv = new_argv

import refetch_orbitrap_hits
refetch_orbitrap_hits.main()
