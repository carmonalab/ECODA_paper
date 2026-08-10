#!/bin/bash
#
# Best-effort sync-status email helper (login-node scripts only).
#
# Source AFTER slurm_config.sh (uses USER_EMAIL, exported there):
#   source "$(dirname "${BASH_SOURCE[0]}")/../utils/bash/sync_status_email.sh"
#
# Provided:
#   notify_sync_status <subject> <body>
#       Sends a best-effort email to ${USER_EMAIL} via the cluster's existing
#       MTA (the same one that already delivers Slurm's --mail-user job
#       emails). Mail binary discovery order:
#         1. NOTIFY_MAIL_BIN env override (test hook; e.g. NOTIFY_MAIL_BIN=cat
#            prints the body to stdout instead of sending)
#         2. first available of mailx, mail, sendmail via `command -v`
#       If no binary is found the email is skipped with a stdout note and the
#       function returns 0 — it never fails the caller. sendmail gets the
#       message (To:/Subject: headers + body) piped to `sendmail -t`;
#       mailx/mail get `-s <subject> ${USER_EMAIL}` with the body on stdin;
#       any other (test) binary gets the body on stdin.
#
#   The email is inherently decoupled from Slurm's own job emails (those fire
#   at job END/FAIL; this one fires at sync time, potentially minutes later
#   after the rsync) — its unique value is the answer to "did the NAS sync
#   actually complete".
#
#   array_task_report <job_id>
#       Prints per-task rows of the array job as TASKID<TAB>STATE<TAB>ELAPSED
#       <TAB>EXITCODE (one line per task, numeric-sorted by task id), for rows
#       matching ^<job_id>_[0-9]+ (via `sacct -j <id> -n -P
#       --format=JobID,State,Elapsed,ExitCode`). The regex excludes the
#       .batch/.extern/.0 step rows and the bare array master row (master
#       Elapsed = array wall time, captured separately via array_wall_time).
#       Prints nothing when sacct is unavailable (command -v guard) or returns
#       no matching rows; callers embed the output in email bodies and map
#       task ids to datasets/methods.
#   array_wall_time <job_id>
#       Prints the array master row's Elapsed (`sacct -j <id> -X -n
#       --format=Elapsed`), else "n/a" (sacct unavailable/purged/unknown id).
#   Both are safe under `set -euo pipefail` (guarded pipelines, always exit 0)
#   and print to stdout for command substitution into email bodies.

notify_sync_status() {
  local SUBJECT="$1"
  local BODY="$2"
  local MAIL_BIN="${NOTIFY_MAIL_BIN:-}"

  if [[ -z "${MAIL_BIN}" ]]; then
    for candidate in mailx mail sendmail; do
      if command -v "${candidate}" >/dev/null 2>&1; then
        MAIL_BIN="${candidate}"
        break
      fi
    done
  fi

  if [[ -z "${MAIL_BIN}" ]]; then
    echo "NOTE: no mail binary (mailx/mail/sendmail) found; skipping sync-status email."
    return 0
  fi

  local STATUS=0
  case "$(basename "${MAIL_BIN}")" in
    sendmail)
      {
        echo "To: ${USER_EMAIL}"
        echo "Subject: ${SUBJECT}"
        echo
        printf '%s\n' "${BODY}"
      } | "${MAIL_BIN}" -t || STATUS=$?
      ;;
    mailx|mail)
      printf '%s\n' "${BODY}" | "${MAIL_BIN}" -s "${SUBJECT}" "${USER_EMAIL}" || STATUS=$?
      ;;
    *)
      printf '%s\n' "${BODY}" | "${MAIL_BIN}" || STATUS=$?
      ;;
  esac

  if [[ ${STATUS} -ne 0 ]]; then
    echo "WARNING: sync-status email via ${MAIL_BIN} failed (exit ${STATUS}); continuing." >&2
    return 0
  fi
  echo "Sync-status email sent to ${USER_EMAIL} via ${MAIL_BIN}."
}

# ---------------------------------------------------------------------------
# Job-duration report helpers (sacct Elapsed, HH:MM:SS as-is; see header)
# ---------------------------------------------------------------------------

array_task_report() {
  local JOB_ID="$1"
  if ! command -v sacct >/dev/null 2>&1; then
    return 0
  fi
  sacct -j "${JOB_ID}" -n -P --format=JobID,State,Elapsed,ExitCode 2>/dev/null | \
  awk -F'|' -v jid="${JOB_ID}" '
    {
      jobid = $1
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", jobid)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4)
      if (jobid ~ "^" jid "_[0-9]+$") {
        printf "%s\t%s\t%s\t%s\n", substr(jobid, length(jid) + 2), $2, $3, $4
      }
    }' | sort -t "$(printf '\t')" -k1,1n || true
}

array_wall_time() {
  local JOB_ID="$1"
  local WALL="n/a"
  if command -v sacct >/dev/null 2>&1; then
    WALL="$(sacct -j "${JOB_ID}" -X -n --format=Elapsed 2>/dev/null | head -1 | tr -d '[:space:]' || true)"
    [[ -n "${WALL}" ]] || WALL="n/a"
  fi
  printf '%s' "${WALL}"
}
