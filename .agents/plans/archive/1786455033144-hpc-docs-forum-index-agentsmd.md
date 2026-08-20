# Index HPC docs & forum for semantic search + AGENTS.md corrections

## Goal

1. Fix the factually wrong "stale NFS views" wording in AGENTS.md — the doc crawl established `$HOME` is **BeeGFS** (see `docs/hpc_docs/storage_on_hpc.md` after this plan runs).
2. Download the UNIGE HPC official wiki (`/hpc/` namespace, ~21 pages per `?do=index`) and the community forum (Discourse, topics from 2024-01-01 + pinned) as Markdown into `docs/hpc_docs/` and `docs/hpc_forum/`, **committed to the repo**, so the built-in `semantic_search` tool can retrieve cluster knowledge locally (live web searches are often dropped by agents for efficiency).
3. Add a "knowledge base retrieval protocol" section to AGENTS.md instructing agents to `semantic_search` the local snapshots before guessing cluster specifics, with live-site URLs as fallback.

## Decisions (confirmed with user)

- **Scope**: all wiki pages under `/hpc/` + forum topics `created_at >= 2024-01-01` + globally pinned topics regardless of date + `forum_index.md` listing **all** topics (id | title | created | posts | pinned | url).
- **Commit to git** (fresh clones and the HPC clone get the snapshot; Kilo indexes tracked workspace files reliably). Expected ~1–2k files, ~5–10 MB.
- **Scraper**: stdlib-only Python (urllib + html.parser + json) committed at `src/utils/py/fetch_hpc_docs.py`, re-runnable. **No pixi.toml/pixi.lock changes** (guardrail).

## Files touched

- `AGENTS.md` (2 edits, below)
- `src/utils/py/fetch_hpc_docs.py` (new)
- `docs/hpc_docs/*.md` + `docs/hpc_docs/README.md` (new, generated + manual)
- `docs/hpc_forum/*.md` + `docs/hpc_forum/forum_index.md` + `docs/hpc_forum/README.md` (new, generated + manual)
- This plan → archived after implementation (repo Task Completion Workflow)

## Tasks (ordered)

### 1. AGENTS.md — BeeGFS wording fix (line 141)

Replace `worker nodes can serve stale NFS views` with
`worker nodes can serve stale BeeGFS client-cache views ($HOME is BeeGFS, not NFS — see docs/hpc_docs/storage_on_hpc.md)`.
Keep the rest of the sentence unchanged.

### 2. AGENTS.md — new section `## HPC knowledge base (offline, indexed)`

Insert between `## HPC general information` and `## Batch effect analysis dataset info`:

```markdown
## HPC knowledge base (offline, indexed)
- Local snapshots of the UNIGE HPC documentation and community forum are committed
  here: `docs/hpc_docs/` (official wiki, /hpc/ namespace) and `docs/hpc_forum/`
  (forum topics from 2024-01-01, plus pinned; full topic list in
  `docs/hpc_forum/forum_index.md`). Snapshot date: 2026-08-11 — content may drift
  from the live sites. Refresh: `python3 src/utils/py/fetch_hpc_docs.py docs` and
  `... forum`, then commit.
- Retrieval protocol: before guessing cluster specifics (partitions, modules,
  quotas, storage semantics, known issues), run a `semantic_search` query targeted
  at `docs/hpc_docs` and `docs/hpc_forum` (e.g. "BeeGFS attribute cache
  staleness", "current issues on HPC cluster", "debug-cpu partition limits").
  Fall back to the live sites only if results are ambiguous:
  https://doc.eresearch.unige.ch/hpc/start and https://hpc-community.unige.ch/
  (2026 "current issues" thread: https://hpc-community.unige.ch/t/4185).
- The forum is Discourse; its JSON API (/latest.json, /t/<id>.json) is used by the
  scraper and returns post markdown in the `raw` field.
```

### 3. Write `src/utils/py/fetch_hpc_docs.py`

- **stdlib only**: `urllib.request`, `html.parser`, `json`, `re`, `argparse`, `pathlib`, `time`, `datetime`. Runs on macOS python3 (no pixi env needed).
- **CLI**: subcommands `docs` and `forum`; global `--out-dir` (default `docs/`) and `--force` (re-write existing files; default: skip unchanged-URL files? No — always overwrite; snapshot refresh must be a full re-scrape).
- **`docs` subcommand**:
  - Seed `https://doc.eresearch.unige.ch/hpc/start`; BFS over same-host links whose path starts with `/hpc/` (strip `#fragments`, skip `?do=` export links).
  - Extract content from `<main>` (fallback: `div#dokuwiki__content`, then `<body>`).
  - Minimal HTML→Markdown converter via `HTMLParser`: h1–h6 → `#…`, `pre/code` → fenced block, `li` → `- `, `tr/td` → pipe rows, `hr` → `---`, paragraphs → blank lines, images → alt text. Fidelity over completeness: semantic_search needs readable text, not perfect Markdown.
  - Write `docs/hpc_docs/<path-slug>.md` (e.g. `storage_on_hpc.md`) with header `# Source: <url>` + `# Snapshot: 2026-08-11`.
  - Print per-page failures, continue; exit 1 if zero pages written; summary count at end.
- **`forum` subcommand**:
  - Paginate `https://hpc-community.unige.ch/latest.json?page=N` (start 0, stop when `topic_list.topics` empty or `more_topics_url` stops advancing; cap at e.g. 400 pages as a safety valve).
  - Write `docs/hpc_forum/forum_index.md`: all topics, `id | title | created | posts | pinned | url`, sorted by id descending.
  - Select topics: `created_at >= 2024-01-01` OR `pinned_globally == true` (or `pinned == true`).
  - For each selected topic: `GET /t/<id>.json`; if `posts_count > 30`, page through additional posts (Discourse topic JSON `?page=N`; loop until all post numbers fetched).
  - Write `docs/hpc_forum/t<id>-<slug>.md`: header (title, source URL, created, tags/category, snapshot date), then per post `## Post <n> by @<username> (<date>)` + the post's `raw` field (already Markdown).
  - Politeness: `User-Agent: ECODA-agent-doc-indexer/1.0`, ~0.5 s sleep between topic fetches, 2 retries with backoff on 429/5xx.
  - Summary: topics selected/fetched/failed, total size.

### 4. Run the scraper locally (macOS)

```bash
python3 src/utils/py/fetch_hpc_docs.py docs
python3 src/utils/py/fetch_hpc_docs.py forum   # ~15–30 min; run once, locally
```

Verify: docs ≈ 21–40 files; forum ≈ 1–2k files; total < ~15 MB (`du -sh docs/hpc_docs docs/hpc_forum`).

### 5. Manual READMEs

- `docs/hpc_docs/README.md` and `docs/hpc_forum/README.md`: source URLs, snapshot date 2026-08-11, crawl scope, re-scrape command, drift warning.

### 6. Validation

1. `grep -n "BeeGFS" AGENTS.md` → new wording + KB section present.
2. Spot-check converted files: `docs/hpc_docs/storage_on_hpc.md` contains "BeeGFS"; a forum thread md contains readable post markdown; no raw `<div>`/`<span>` HTML leakage.
3. **Semantic search probe** (the whole point): in this workspace run `semantic_search` with e.g. "BeeGFS storage home directory" and "current issues HPC cluster 2026"; confirm hits under `docs/hpc_docs` / `docs/hpc_forum`. If Kilo does not index the files, document the fallback (grep + live URLs) in AGENTS.md — but do not leave the section claiming functionality that does not work.
4. `git status` shows only: AGENTS.md, the new script, `docs/hpc_docs/**`, `docs/hpc_forum/**`. **No pixi.toml / pixi.lock changes.**

### 7. Commit & push (repo Task Completion Workflow)

Archive this plan to `.kilo/plans/archive/`, `git add .`, commit (e.g. `Index HPC docs + forum for semantic search; AGENTS.md BeeGFS correction + KB protocol`), push.

## Failure modes & notes

- **Forum rate-limiting/blocking**: retries + sleep; if persistent, reduce scope via `--since 2025-01-01` and note in the forum README.
- **Crawl duration**: ~1–2k topics × ~0.5–1 s ≈ 15–30 min; run in a single local session (no HPC involvement; network is fine — both sites were reachable from this machine during planning).
- **Repo bloat**: if `du -sh` exceeds ~15 MB, stop and re-scope (drop pre-2025 topics) before committing.
- **DokuWiki layout drift**: content extraction falls back body→generic; degraded but acceptable HTML→text is fine for search.
- **Out of scope**: NASAC/other wiki namespaces, `.kilocode/rules` file (user explicitly dropped it), ARCHITECTURE.md/README.md updates.
