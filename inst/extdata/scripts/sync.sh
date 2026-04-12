#!/bin/bash

SRC=${1:-~/workspace/project}
DST=${2:-user@host:/remote/project}
PATTERN=${3:-'r.*.r'}

LOCKFILE=/tmp/rsync-sync.lock
PIDFILE=/tmp/rsync-sync.pid

# Expand ~ in path
SRC=$(eval echo "$SRC")

if [ ! -d "$SRC" ]; then
  notify-send "SRC directory does not exist: $SRC"
  exit 1
fi

# ========= Startup check =========

if [ -f "$PIDFILE" ]; then
  OLD_PID=$(cat "$PIDFILE")

  if ps -p "$OLD_PID" > /dev/null 2>&1; then
    notify-send "Sync is already running (PID=$OLD_PID), exiting"
    exit 1
  else
    notify-send "Stale PID file found, cleaning up"
    rm -f "$PIDFILE"
  fi
fi

echo $$ > "$PIDFILE"

# Cleanup on exit
cleanup() {
  rm -f "$LOCKFILE" "$PIDFILE"
  exit
}
trap cleanup INT TERM EXIT

notify-send "Sync started: $SRC -> $DST"

# ========= Main loop =========

while true; do
  inotifywait -r -e modify,create,delete,move . >/dev/null 2>&1

  # Debounce to avoid excessive triggering
  sleep 1

  # Prevent concurrent rsync execution
  if [ -f "$LOCKFILE" ]; then
    continue
  fi

  touch "$LOCKFILE"

  rsync -az \
    --include='*/' \
    --include="$PATTERN" \
    --exclude='*' \
    "$SRC"/ "$DST"

  rm -f "$LOCKFILE"
done
