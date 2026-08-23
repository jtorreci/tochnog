#!/bin/bash
# Lanza opencode bajo vigilancia para auditar salidas bruscas.
#
# Sustituye a `opencode` (alias: alias opencode=/ruta/scripts/opencode-watch.sh).
# Mientras corre, un latido (heartbeat) se apunta en session-audit.log cada 60s
# con uptime + memoria. Si WSL/opencode se cae (suspension del portatil, OOM...),
# la ultima linea del log queda sin el cierre "termino"; al arrancar de nuevo,
# el script detecta la salida brusca y lo deja anotado en el propio log.
#
# Uso: ./scripts/opencode-watch.sh [args-de-opencode...]

set -uo pipefail

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOG="$REPO_DIR/session-audit.log"
INTERVAL=60

log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')" "$*" >> "$LOG"; }

# --- Detectar si la sesion anterior murio sin cerrarse limpio --------------
LAST="$(tail -n 1 "$LOG" 2>/dev/null || true)"
if printf '%s' "$LAST" | grep -q 'WATCH.*HB alive'; then
  log "WATCH: !! SESION ANTERIOR ABRUPTA (ultima linea = latido sin cierre limpio)"
  log "WATCH:     ultima actividad antes del corte: $LAST"
fi

log "WATCH: arrancando opencode (branch $(git -C "$REPO_DIR" branch --show-current 2>/dev/null || echo '?'), uptime $(uptime -p))"

# --- Heartbeat en segundo plano --------------------------------------------
(
  while true; do
    FREE="$(awk '/MemFree|MemAvailable/ {printf "%s=%dKB ", $1, $2}' /proc/meminfo)"
    log "WATCH: HB alive uptime=$(cut -d. -f1 /proc/uptime)s load=$(cut -d' ' -f1 /proc/loadavg) $FREE"
    sleep "$INTERVAL"
  done
) &
HB_PID=$!
trap 'kill "$HB_PID" 2>/dev/null || true' EXIT

opencode "$@"
CODE=$?
log "WATCH: opencode termino con codigo $CODE"
exit "$CODE"
