#!/usr/bin/env python3
"""Fetch UNIGE HPC documentation (DokuWiki /hpc/ namespace) and community forum
(Discourse) into local Markdown snapshots under docs/, so agents can retrieve
cluster knowledge via semantic_search without hitting the live sites.

Stdlib only (urllib, html.parser, json, re, argparse, pathlib, time) — runs on
any python3 (macOS system python3); pixi.toml/pixi.lock must stay untouched.

Usage:
  python3 src/utils/py/fetch_hpc_docs.py docs   [--out-dir docs]
  python3 src/utils/py/fetch_hpc_docs.py forum  [--out-dir docs] [--since 2024-01-01] [--max-pages 400]

Re-scraping always overwrites existing files (full snapshot refresh).
"""

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

USER_AGENT = "ECODA-agent-doc-indexer/1.0"
SNAPSHOT_DATE = "2026-08-11"
DOCS_BASE = "https://doc.eresearch.unige.ch"
FORUM_BASE = "https://hpc-community.unige.ch"
DEFAULT_SINCE = "2024-01-01"
DEFAULT_MAX_PAGES = 400
SLEEP_TOPIC = 0.5
RETRIES = 2


def http_get(url, timeout=90):
    """GET url with User-Agent; 2 retries with backoff on 429/5xx/timeouts."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < RETRIES:
                time.sleep(5 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError:
            if attempt < RETRIES:
                time.sleep(5 * (attempt + 1))
                continue
            raise


def iso_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_md(path, header_lines, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    content = header_lines + ["---", ""] + body
    path.write_text("\n".join(content).rstrip() + "\n", encoding="utf-8")


def safe_slug(name):
    return re.sub(r"[^a-zA-Z0-9._-]+", "_", name.strip()).strip("_")


# ---------------------------------------------------------------------------
# DokuWiki raw markup -> Markdown (docs subcommand)
# ---------------------------------------------------------------------------

DOKU_HEADING = re.compile(r"^(={2,7})\s*(.*?)\s*\1\s*$")
DOKU_LINK = re.compile(r"\[\[([^\]|]+)(?:\|([^\]]*))?\]\]")
DOKU_MONO = re.compile(r"''([^']+)''")
DOKU_ITALIC = re.compile(r"//([^/\n]+?)//")
DOKU_ITALIC_LINK = re.compile(r"//(\[[^\n]*?\]\([^\n]*?\))//")
DOKU_IMAGE = re.compile(r"\{\{[^}]*\}\}")
DOKU_TAG = re.compile(r"</?(?:sub|sup|del|ins|u|nowiki|span|div|note|wrap|WRAP|html)[^>]*>")


def _doku_inline(line):
    """Inline DokuWiki -> Markdown transforms for a single text line."""
    line = DOKU_LINK.sub(lambda m: _doku_link(m), line)
    line = DOKU_MONO.sub(r"`\1`", line)
    line = DOKU_ITALIC_LINK.sub(r"*\1*", line)
    line = DOKU_ITALIC.sub(r"*\1*", line)
    line = DOKU_TAG.sub("", line)
    line = re.sub(r"\\\\", " ", line)
    line = re.sub(r"\(\((?:[^()]|\([^)]*\))*\)\)", lambda m: m.group(0)[2:-2] or "", line)
    return line.strip()


def dokuwiki_to_markdown(text):
    """Line-based DokuWiki -> Markdown pass. Fidelity over completeness."""
    lines = text.split("\n")
    out = []
    in_code = None
    for raw in lines:
        line = raw.rstrip()
        cm = re.match(r"^\s*<code(?:\s+([^>]*))?>\s*(.*)$", line)
        if in_code is not None:
            if "</code>" in line:
                code_body = line.split("</code>", 1)[0]
                out.append(code_body)
                out.append("```")
                in_code = None
            else:
                out.append(line)
            continue
        if cm:
            lang = cm.group(1).strip() if cm.group(1) else ""
            out.append("```" + (lang.split()[0] if lang else ""))
            if cm.group(2):
                out.append(cm.group(2))
            in_code = True
            continue
        if line.lstrip().startswith("<file"):
            out.append("```")
            if "</file>" in line:
                out.append(line.split("</file>", 1)[0].split(">", 1)[1])
                out.append("```")
            else:
                in_code = True
            continue
        if in_code:
            out.append(line)
            if "</file>" in line:
                out.append("```")
                in_code = False
            continue
        if re.match(r"^\s*<html>", line):
            in_code = True
            continue
        stripped = line.strip()
        if not stripped:
            out.append("")
            continue
        if stripped == "\\":
            continue
        if stripped.startswith("+++"):
            continue
        if re.match(r"^[?!]{3}\s", stripped):
            prefix = "**Q:** " if stripped.startswith("?") else "**A:** "
            out.append(prefix + _doku_inline(stripped[3:].strip()))
            continue
        if stripped.startswith("{{") and stripped.endswith("}}"):
            if stripped.startswith("{{METATOC"):
                continue
            m = DOKU_IMAGE.match(stripped)
            if m:
                continue
        if re.match(r"^-{4,}$", stripped):
            out.append("---")
            continue
        note = re.match(r"^<note(?: [^>]*)?>(.*)$", stripped)
        if note:
            out.append("> " + _doku_inline(note.group(1)))
            continue
        if stripped.startswith("</note>"):
            continue
        if re.match(r"^=(\?)?=+", stripped) and re.search(r"=+(\?)?=$", stripped):
            stripped = re.sub(r"^=(\?)?", "=", stripped)
            stripped = re.sub(r"(\?)?=$", "=", stripped)
            stripped = stripped.strip()
        m = DOKU_HEADING.match(stripped)
        if m:
            level = 8 - len(m.group(1))
            out.append("#" * level + " " + m.group(2).strip())
            continue
        if stripped.startswith(("^", "|")) and not stripped.startswith(("^^", "||")):
            cells = [c.strip() for c in stripped.split("|")[1:-1]] \
                if stripped.startswith("|") else \
                [c.strip() for c in stripped.split("^")[1:-1]]
            if any(c.startswith(":::") for c in cells):
                prev = out[-1] if out and out[-1].startswith("|") else "| |"
                pcell = prev.split("|")[1:-1]
                cells = [pcell[i] if (c.startswith(":::") and i < len(pcell)) else c
                         for i, c in enumerate(cells)]
            cells = [_doku_inline(c) for c in cells]
            if stripped.startswith("^"):
                out.append("| " + " | ".join(cells) + " |")
                out.append("| " + " | ".join(["---"] * len(cells)) + " |")
            else:
                out.append("| " + " | ".join(cells) + " |")
            continue
        lm = re.match(r"^(\s{2,})([*\-])\s+(.*)$", line)
        if lm:
            indent = "  " * max(0, len(lm.group(1)) // 2 - 1)
            out.append(indent + "- " + _doku_inline(lm.group(3)))
            continue
        out.append(_doku_inline(line.rstrip("\\")))
    body = "\n".join(out)
    body = re.sub(r"\n{3,}", "\n\n", body)
    return body.strip() + "\n"


def _doku_link(m):
    target, label = m.group(1).strip(), (m.group(2) or "").strip()
    if "://" in target:
        url = target
    else:
        page, _, anchor = target.partition("#")
        path = page.replace(":", "/")
        url = f"{DOCS_BASE}/{path}"
        if anchor:
            url += f"#{anchor}"
    if not label:
        label = url
    return f"[{label}]({url})"


def docs_crawl(out_dir):
    out_dir = Path(out_dir)
    pages = _docs_seed_pages()
    if not pages:
        print("ERROR: no pages found in sitemap", file=sys.stderr)
        return 1
    queue = sorted(pages)
    seen = {p.lower() for p in pages}
    written = 0
    failures = []
    while queue:
        page = queue.pop(0)
        url = f"{DOCS_BASE}/_export/raw/{urllib.parse.quote(page, safe=':')}"
        try:
            markup = http_get(url).decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 403:
                print(f"  SKIP {page}: HTTP 403 (alias/denied page)", file=sys.stderr)
                continue
            failures.append((page, f"{type(e).__name__}: {e}"))
            print(f"  FAIL {page}: {e}", file=sys.stderr)
            continue
        except Exception as e:
            failures.append((page, f"{type(e).__name__}: {e}"))
            print(f"  FAIL {page}: {e}", file=sys.stderr)
            continue
        for m in re.findall(r"\[\[(hpc|wiki|linux):([^\]|#]+)", markup):
            found = re.sub(r":+", ":", f"{m[0]}:{m[1]}").rstrip(":")
            if not found.startswith("hpc:") or not found[4:]:
                continue
            if found.lower() not in seen:
                seen.add(found.lower())
                queue.append(found)
        body = dokuwiki_to_markdown(markup)
        slug = page.split(":", 1)[1] if ":" in page else page
        rel = Path(*slug.split("/")) if "/" in slug else Path(slug)
        path = out_dir / "hpc_docs" / rel.with_suffix(".md")
        write_md(
            path,
            [f"# Source: {DOCS_BASE}/{page.replace(':', '/')}", f"# Snapshot: {SNAPSHOT_DATE}", f"# Crawled: {iso_now()}", ""],
            body.split("\n"),
        )
        written += 1
        time.sleep(0.3)
    print(f"docs: {written} pages written, {len(failures)} failed")
    for page, err in failures:
        print(f"  - {page}: {err}", file=sys.stderr)
    return 0 if written else 1


def _docs_seed_pages():
    try:
        html = http_get(f"{DOCS_BASE}/hpc/start?do=index").decode("utf-8", errors="replace")
        pages = sorted({p for p in re.findall(r'href="/hpc/([^"?#]+)"', html)})
    except Exception as e:
        print(f"WARNING: sitemap fetch failed ({e}); seeding from /hpc/start only", file=sys.stderr)
        pages = ["hpc:start"]
    if "hpc:start" not in pages:
        pages.append("hpc:start")
    return pages


# ---------------------------------------------------------------------------
# Discourse -> Markdown (forum subcommand)
# ---------------------------------------------------------------------------

class DiscourseHTMLToMD(HTMLParser):
    """Convert Discourse `cooked` HTML into readable Markdown."""

    BLOCK_END = {"p", "div", "h1", "h2", "h3", "h4", "h5", "h6", "li",
                 "blockquote", "pre", "tr", "ul", "ol", "table", "section",
                 "aside", "article", "hr", "br", "table", "thead", "tbody"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.lines = []
        self.buf = ""
        self.list_depth = 0
        self.quote_depth = 0
        self.pre = False
        self.pre_buf = []
        self.heading = 0
        self.link_href = None
        self.link_text = []
        self.img = []

    def _flush(self):
        text = self.buf.strip()
        if not text:
            self.buf = ""
            return
        prefix = "> " * self.quote_depth
        if self.heading:
            text = "#" * self.heading + " " + text
        elif self.list_depth:
            text = "  " * (self.list_depth - 1) + "- " + text
        self.lines.append(prefix + text)
        self.buf = ""

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in ("p", "div", "li", "tr", "ul", "ol", "table", "section", "aside", "article", "thead", "tbody"):
            self._flush()
            if tag in ("ul", "ol"):
                self.list_depth += 1
        elif tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self._flush()
            self.heading = int(tag[1])
        elif tag == "blockquote":
            self._flush()
            self.quote_depth += 1
        elif tag == "pre":
            self._flush()
            self.pre = True
        elif tag == "code":
            if not self.pre:
                self.buf += "`"
        elif tag == "br":
            self._flush()
        elif tag == "hr":
            self._flush()
            self.lines.append("---")
        elif tag == "a":
            href = a.get("href", "")
            is_anchor = bool(re.match(r"^#(p-)?\d+", href)) or "name=" in href or "aria-label" in a
            if href and not is_anchor:
                if href.startswith("/"):
                    href = f"{FORUM_BASE}{href}"
                self.link_href = href
                self.link_text = []
        elif tag == "img":
            alt = a.get("alt", "") or a.get("title", "")
            if alt:
                self.buf += alt

    def handle_endtag(self, tag):
        if tag in ("p", "div", "li", "tr", "table", "section", "aside", "article", "thead", "tbody"):
            self._flush()
        elif tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self._flush()
            self.heading = 0
        elif tag == "blockquote":
            self._flush()
            self.quote_depth = max(0, self.quote_depth - 1)
        elif tag in ("ul", "ol"):
            self._flush()
            self.list_depth = max(0, self.list_depth - 1)
        elif tag == "pre":
            self._flush()
            self.pre = False
            code = "".join(self.pre_buf).strip("\n")
            self.pre_buf = []
            if code:
                self.lines.append("```")
                self.lines.extend(code.split("\n"))
                self.lines.append("```")
        elif tag == "code":
            if not self.pre:
                self.buf += "`"
        elif tag == "a":
            if self.link_href is not None:
                text = "".join(self.link_text).strip()
                url = self.link_href
                self.link_href = None
                self.link_text = []
                if text:
                    self.buf += f"[{text}]({url})"

    def handle_data(self, data):
        if self.pre:
            self.pre_buf.append(data)
            return
        if self.link_href is not None:
            self.link_text.append(data)
        self.buf += data


def cooked_to_markdown(cooked):
    p = DiscourseHTMLToMD()
    p.feed(cooked)
    p.close()
    p._flush()
    body = "\n".join(p.lines)
    body = re.sub(r"\n{3,}", "\n\n", body)
    return body.strip() + "\n"


def forum_posts(topic_id, posts_count, max_pages):
    """Fetch all posts of a topic (paginated). Returns list of post dicts."""
    posts = {}
    page = 0
    while page <= max_pages:
        url = f"{FORUM_BASE}/t/{topic_id}.json"
        if page:
            url += f"?page={page}"
        try:
            raw = http_get(url)
        except Exception as e:
            print(f"  FAIL topic {topic_id} page {page}: {e}", file=sys.stderr)
            break
        data = json.loads(raw.decode("utf-8", errors="replace"))
        stream = data.get("post_stream", {}).get("posts", [])
        new = False
        for p in stream:
            if p.get("post_number") not in posts:
                posts[p.get("post_number")] = p
                new = True
        if not stream or (posts_count and len(posts) >= posts_count) or not new:
            break
        page += 1
        time.sleep(SLEEP_TOPIC)
    return [posts[n] for n in sorted(posts)]


def forum_crawl(out_dir, since, max_pages):
    out_dir = Path(out_dir)
    since_dt = datetime.strptime(since, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    topics = []
    page = 0
    while page <= max_pages:
        url = f"{FORUM_BASE}/latest.json?page={page}"
        try:
            data = json.loads(http_get(url).decode("utf-8", errors="replace"))
        except Exception as e:
            print(f"FAIL latest.json page {page}: {e}", file=sys.stderr)
            page += 1
            continue
        batch = data.get("topic_list", {}).get("topics", [])
        if not batch:
            break
        seen_ids = {t["id"] for t in topics}
        topics.extend(t for t in batch if t["id"] not in seen_ids)
        page += 1
        time.sleep(SLEEP_TOPIC)
    topics.sort(key=lambda t: t["id"], reverse=True)
    print(f"forum: {len(topics)} topics listed (pages {page})")

    index_lines = ["# HPC community forum — topic index",
                   "",
                   f"# Snapshot: {SNAPSHOT_DATE}",
                   "",
                   "All topics visible in /latest.json, sorted by id descending.",
                   "",
                   "| id | title | created | posts | pinned | url |",
                   "| --- | --- | --- | --- | --- | --- |"]
    selected = []
    for t in topics:
        created = (t.get("created_at") or "")[:10]
        pinned = bool(t.get("pinned_globally") or t.get("pinned"))
        title = t.get("title", "?").replace("|", "\\|")
        index_lines.append(
            f"| {t['id']} | {title} | {created} | {t.get('posts_count', '?')} | "
            f"{'yes' if pinned else ''} | {FORUM_BASE}/t/{t['id']} |")
        if created >= since or pinned:
            selected.append(t)
    index_path = out_dir / "hpc_forum" / "forum_index.md"
    index_path.parent.mkdir(parents=True, exist_ok=True)
    index_path.write_text("\n".join(index_lines) + "\n", encoding="utf-8")
    print(f"forum: {len(selected)} topics selected (created >= {since} or pinned)")

    written = 0
    failed = []
    for t in selected:
        tid = t["id"]
        slug = safe_slug(t.get("slug", "")) or str(tid)
        fname = f"t{tid}-{slug}.md"
        path = out_dir / "hpc_forum" / fname
        try:
            posts = forum_posts(tid, t.get("posts_count", 0), 10)
            if not posts:
                failed.append((tid, "no posts"))
                continue
            header = [f"# {t.get('title', '?')}", "",
                      f"- Source: {FORUM_BASE}/t/{tid}", "",
                      f"- Created: {t.get('created_at', '?')}", ""]
            if t.get("tags"):
                tag_names = [x.get("name", str(x)) if isinstance(x, dict) else str(x)
                             for x in t["tags"]]
                header += [f"- Tags: {', '.join(tag_names)}", ""]
            header += [f"- Posts: {t.get('posts_count', len(posts))}", "",
                       f"- Category: {t.get('category_id', '?')}", "",
                       f"- Pinned: {bool(t.get('pinned_globally') or t.get('pinned'))}", "",
                       f"- Snapshot: {SNAPSHOT_DATE}", ""]
            body = []
            for p in posts:
                body += [f"## Post {p.get('post_number', '?')} by @{p.get('username', '?')} "
                         f"({p.get('created_at', '?')})", ""]
                body += cooked_to_markdown(p.get("cooked", "")).split("\n") + [""]
            write_md(path, header, body)
            written += 1
        except Exception as e:
            failed.append((tid, f"{type(e).__name__}: {e}"))
            print(f"  FAIL t{tid}: {e}", file=sys.stderr)
        time.sleep(SLEEP_TOPIC)
    print(f"forum: {written} topics written, {len(failed)} failed")
    for tid, err in failed:
        print(f"  - t{tid}: {err}", file=sys.stderr)
    return 0 if written else 1


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("subcommand", choices=["docs", "forum"])
    ap.add_argument("--out-dir", default="docs", help="output root (default: docs)")
    ap.add_argument("--force", action="store_true",
                    help="accepted for explicitness; re-scrapes always overwrite")
    ap.add_argument("--since", default=DEFAULT_SINCE, help=f"forum: only topics created >= this date (default: {DEFAULT_SINCE})")
    ap.add_argument("--max-pages", type=int, default=DEFAULT_MAX_PAGES,
                    help=f"forum: safety cap on list pages (default: {DEFAULT_MAX_PAGES})")
    args = ap.parse_args(argv)
    t0 = time.time()
    if args.subcommand == "docs":
        rc = docs_crawl(args.out_dir)
    else:
        rc = forum_crawl(args.out_dir, args.since, args.max_pages)
    print(f"done in {time.time() - t0:.0f}s")
    return rc


if __name__ == "__main__":
    sys.exit(main())
