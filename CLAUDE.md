@AGENTS.md

## Claude Code specifics

- The skills above are discovered from `.claude/skills/` in this repo; invoke one
  explicitly with `/mutadock-setup`, `/mutadock-mutation`, and so on.
- `.claude/skills/` is a **generated mirror** of `.agents/skills/` (which Codex
  reads). Edit `.agents/skills/`, then run `python scripts/sync_agent_skills.py`.
  `pytest tests/test_skills` fails if the two trees drift apart.
- ΔΔG and docking commands are long-running and CPU-bound. Run them in the
  background or with a pilot-sized workload; never start a default-size
  `md_mutate` without saying so first.
