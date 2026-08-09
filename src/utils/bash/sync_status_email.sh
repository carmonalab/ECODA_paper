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
