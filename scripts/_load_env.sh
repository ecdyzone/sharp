#!/usr/bin/env bash
# Load the project `.env` into the environment. Source me; don't execute me.
#
#   source "$(dirname "${BASH_SOURCE[0]}")/_load_env.sh"
#
# `set -a` auto-exports every assignment, so values reach child processes that
# read the environment themselves (e.g. `deepbgc download` and
# DEEPBGC_DOWNLOADS_DIR) without the caller plumbing them through explicitly.
#
# Sourcing with bash — rather than parsing the file — is what makes the
# `${DATABASES:-$HOME/.local/share}` fallbacks in `.env` work: on the server
# DATABASES comes from the surrounding environment, on a laptop it falls back.
# `src/sharp/config.py` has its own minimal parser for the SHARP_* keys; it does
# not understand `:-`, so keys written that way are for these scripts only.
#
# Missing `.env` is not an error — every key is optional and the callers apply
# their own `: "${VAR:=default}"` fallbacks.

_env_file="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/.env"
if [[ -f "$_env_file" ]]; then
    set -a
    # shellcheck disable=SC1090  # runtime path, nothing to check statically
    source "$_env_file"
    set +a
fi
unset _env_file
