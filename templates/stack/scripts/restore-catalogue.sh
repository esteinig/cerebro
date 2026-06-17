#!/usr/bin/env bash
#
# restore-catalogue.sh (S4-2) — restore a Cerebro catalogue/audit backup.
#
# Operator-gated and destructive: a real restore overwrites live catalogue data,
# so it asks for confirmation (D6). It verifies the archive against the manifest
# BLAKE3 checksum before touching MongoDB, and prints the audit-chain state so you
# can re-verify the chain after loading (a backup can't silently launder a broken
# chain — the manifest records whether the chain was intact at backup time).
#
# Default backend is the filesystem object store (a mounted backup volume / NFS),
# matching the worker's FilesystemObjectStore. The S3 backend (S4-4) will add a
# download step ahead of the same verify/restore flow.
#
# Requirements: bash, jq, b3sum (BLAKE3), and mongorestore (mongodb-database-tools).
#
# Usage:
#   restore-catalogue.sh --uri <mongo-uri> [--backup <id|latest>] [--store <path>]
#                        [--prefix <prefix>] [--dry-run] [--drop]
#
#   --uri      Target MongoDB URI to restore into (required for restore/dry-run).
#   --backup   Backup id (e.g. 20260617T031500Z) or "latest" (default: latest).
#   --store    Object-store root (default: $CEREBRO_BACKUP_STORE_PATH or /backups).
#   --prefix   Key prefix (default: $CEREBRO_BACKUP_PREFIX or "catalogue").
#   --dry-run  Verify checksum + that the archive is loadable (mongorestore
#              --dryRun); never writes.
#   --drop     Drop each collection before restoring it (clean restore).
#
set -euo pipefail

STORE="${CEREBRO_BACKUP_STORE_PATH:-/backups}"
PREFIX="${CEREBRO_BACKUP_PREFIX:-catalogue}"
BACKUP_ID="latest"
TARGET_URI=""
DRY_RUN=0
DROP=0

die() { echo "error: $*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --uri)     TARGET_URI="${2:-}"; shift 2 ;;
    --backup)  BACKUP_ID="${2:-}"; shift 2 ;;
    --store)   STORE="${2:-}"; shift 2 ;;
    --prefix)  PREFIX="${2:-}"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --drop)    DROP=1; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

command -v jq          >/dev/null || die "jq is required"
command -v b3sum        >/dev/null || die "b3sum (BLAKE3) is required"
command -v mongorestore >/dev/null || die "mongorestore (mongodb-database-tools) is required"
[ -n "$TARGET_URI" ] || die "--uri <mongo-uri> is required"

# Resolve "latest": backup ids are sortable timestamps, so the last directory wins.
if [ "$BACKUP_ID" = "latest" ]; then
  BACKUP_ID="$(ls -1 "$STORE/$PREFIX" 2>/dev/null | sort | tail -n1 || true)"
  [ -n "$BACKUP_ID" ] || die "no backups found under $STORE/$PREFIX"
  echo "resolved latest backup: $BACKUP_ID"
fi

DIR="$STORE/$PREFIX/$BACKUP_ID"
ARCHIVE="$DIR/catalogue.archive.gz"
MANIFEST="$DIR/manifest.json"
[ -f "$ARCHIVE" ]  || die "archive not found: $ARCHIVE"
[ -f "$MANIFEST" ] || die "manifest not found: $MANIFEST"

# Verify the archive against the manifest checksum BEFORE any restore.
EXPECTED="$(jq -r '.archive_blake3' "$MANIFEST")"
ACTUAL="$(b3sum --no-names "$ARCHIVE" | awk '{print $1}')"
[ "$EXPECTED" = "$ACTUAL" ] \
  || die "checksum mismatch: manifest=$EXPECTED actual=$ACTUAL (archive corrupt or tampered)"
echo "checksum OK ($ACTUAL)"

CREATED="$(jq -r '.created_at' "$MANIFEST")"
DBLABEL="$(jq -r '.mongo_db' "$MANIFEST")"
CHAIN="$(jq -r '.audit_chain_verified' "$MANIFEST")"
EVENTS="$(jq -r '.audit_event_count' "$MANIFEST")"
echo "backup:    $BACKUP_ID"
echo "created:   $CREATED"
echo "scope:     $DBLABEL"
echo "audit chain verified at backup time: $CHAIN ($EVENTS events)"
if [ "$CHAIN" != "true" ]; then
  echo "WARNING: the audit chain was NOT confirmed intact when this backup was taken." >&2
  echo "         Restoring it will reproduce that state — investigate before relying on it." >&2
fi

if [ "$DRY_RUN" = 1 ]; then
  echo "dry-run: validating the archive is loadable (no writes)…"
  mongorestore --uri="$TARGET_URI" --gzip --archive="$ARCHIVE" --dryRun
  echo "dry-run OK — archive is loadable."
  exit 0
fi

# Destructive restore — operator confirmation gate (D6).
echo
echo "About to restore backup $BACKUP_ID into:"
echo "    $TARGET_URI"
[ "$DROP" = 1 ] && echo "    (collections will be DROPPED before restore)"
echo "This OVERWRITES live catalogue data and cannot be undone."
read -r -p "Type 'restore' to proceed: " CONFIRM
[ "$CONFIRM" = "restore" ] || die "aborted by operator"

if [ "$DROP" = 1 ]; then
  mongorestore --uri="$TARGET_URI" --gzip --archive="$ARCHIVE" --drop
else
  mongorestore --uri="$TARGET_URI" --gzip --archive="$ARCHIVE"
fi

echo
echo "Restore complete. Now re-verify the audit chain (chain-verify-after):"
echo "    GET <api>/files/audit?team=<team>  ->  data.verified must be true"
echo "Compare against the manifest value above (audit_chain_verified=$CHAIN)."
echo "If the live chain now verifies but the manifest said false (or vice-versa),"
echo "stop and investigate before resuming operations."
