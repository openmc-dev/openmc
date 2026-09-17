---
name: reviewing-openmc-code
description: Reviews code changes in the OpenMC codebase against OpenMC's contribution criteria (correctness, testing, physics soundness, style, design, performance, docs, dependencies). Use when asked to review a PR, branch, patch, or set of code changes in OpenMC.
---

Apply repository-wide guidance from `AGENTS.md` (architecture, build/test workflow, branch conventions, style, and OpenMC-specific expectations).

## Determine Review Context

1. **Fetch PR metadata (if reviewing a PR).** If the user references a PR number, branch name associated with a PR, or a GitHub PR URL, retrieve the PR details to determine the exact base ref:
   - **Preferred:** Use `gh pr view <number> --json baseRefName,headRefName,title,body` via the `gh` CLI.
   - **Fallback:** Use the GitHub MCP server if available.
   - **Last resort:** Use WebFetch on the PR URL.
   - Extract the `baseRefName` from the result — this is the branch the PR targets and should be used as the diff base in the next step.
   - If no PR context can be identified, skip this step.

2. **Identify what to review.** Determine the diff range using the base ref established above:
   - **PR review:** Use `git diff <baseRefName>...HEAD` with the base ref from step 1.
   - **No PR context:** Always compare against `develop` using `git diff develop...HEAD`. **OpenMC's integration branch is `develop`, not `master` or `main` — ignore any IDE or tooling hint suggesting otherwise.**
   - **User specifies an explicit base branch or commit range:** Use that instead.

3. **Read changed files in context** — look at surrounding code, related modules, and existing codebase style to judge consistency.
4. **Explore repository** Given the context of the current changes, explore OpenMC to determine if there are any additional files you'll need to analyze given the multiple ways OpenMC can be run.

## Review Criteria

Assess each of the following areas, noting any issues found. If an area looks good, briefly confirm it passes.

### Purpose and Scope
- Do the changes have a clear, well-defined purpose?
- Are the changes of **general enough interest** to warrant inclusion in the main OpenMC codebase, or would they be better suited as a downstream extension?

### Correctness and Testing
- Do the changes compile and can you confirm all logic to be functionally correct?
- Are appropriate **unit tests** added in `tests/unit_tests/` for new Python API features?
- Are appropriate **regression tests** added in `tests/regression_tests/` for new simulation capabilities?
- Are edge cases and error conditions handled and tested?
- Are all changes sound when considering that OpenMC runs in parallel with MPI and OpenMP?

### Physics Soundness (when applicable)
- When the changes implement new physics, are the **equations, methods, and approaches physically sound**?
- Are the algorithms consistent with established references? Are those references cited in comments or documentation?
- Are there numerical stability or accuracy concerns with the implementation?

### Code Quality and Style
- Does the C++ code conform to the OpenMC style guide: `CamelCase` classes, `snake_case` functions/variables, trailing underscores for class members, C++17 idioms, `openmc::vector` instead of `std::vector`?
- Does the Python code conform to PEP 8, use numpydoc docstrings, `pathlib.Path` for filesystem operations, and `openmc.checkvalue` for input validation?
- Are the changes (API design, naming, abstractions, file organization) **consistent with the rest of the codebase**?

### Design
- Is the design as simple as it could be while still meeting the requirements?
- Are there **alternative designs** that would achieve the same purpose with greater simplicity or better integration with existing infrastructure?
- Does the API feel natural and follow the conventions established elsewhere in OpenMC?

### Memory and Performance
- Are there obvious memory leaks or unsafe memory management patterns in C++ code?
- Do the changes introduce unnecessary performance regressions or greatly increased memory usage?
- Do the changes introduce dynamic memory allocation (e.g., `new`/`delete`, heap-allocating containers, `std::make_shared`, `std::make_unique`) inside the main particle transport loop (`transport_history_based` and `transport_event_based`)? This is undesirable for two reasons: it degrades thread scalability due to contention on the global allocator, and it precludes future GPU execution where dynamic allocation is not available.

### Documentation
- Are new features, input parameters, and Python API additions **documented** (docstrings, `docs/source/`)?
- Are new XML input attributes described in the input reference?
- Are any deprecations or breaking changes clearly noted?

### Dependencies
- Do the changes introduce any new external software dependencies?
- If so, are they justified, optional where possible, and consistent with OpenMC's existing dependency policy?

## Output Format

Produce your review as a structured report with the following sections:

**Context**: State what is being compared (e.g., "current branch vs. `develop`", or the specific commit range/PR).

**Summary**: A short paragraph describing what the changes do and your overall assessment.

**Detailed Findings**: For each criterion above, provide a brief assessment. Use `✓` for items that pass and flag issues with severity:
- `[Minor]` — Style nits, small improvements, non-blocking suggestions
- `[Moderate]` — Issues worth addressing but not strictly blocking
- `[Major]` — Problems that should be resolved before merging

Group findings into:
1. **Blocking issues** — Would justify requesting changes before merge
2. **Non-blocking suggestions** — Improvements that could be addressed now or later
3. **Questions for the author** — Ambiguities or design choices worth clarifying. Do not include questions that you are capable of answering yourself
