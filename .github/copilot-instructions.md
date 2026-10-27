# Pull Request Description Guidelines

When drafting or summarizing pull request descriptions, always adhere strictly to the following structural and stylistic rules:

---

### 1. Opening Overview Paragraph
- Open with a clear, synthetic narrative paragraph starting with:
  `This pull request introduces [brief summary of major theme/goal]...`
- State the 2–4 primary areas or systems enhanced (e.g., test runner, CI workflow, documentation).
- Explain *why* the changes were made and their high-level impact.
- Avoid low-level code details or file names in this opening paragraph.

---

### 2. Thematic Grouping (Bold Category Headers)
Group the changes by logical functional domain rather than listing files sequentially. Use standalone bold text with a colon as section headers (e.g., `Smoke test script improvements:`, `Continuous integration (CI) workflow enhancements:`, `Documentation updates:`).

Common category examples (adapt dynamically based on the diff):
- **[Component/Feature] improvements:**
- **Continuous integration (CI) workflow enhancements:**
- **Configuration & dependency updates:**
- **Documentation updates:**

---

### 3. Bullet Point Construction & File Linking
Under each category header, provide 2–4 detailed bullet points conforming to these rules:
- **Tone:** Use active, present-tense phrasing explaining what the system now does (e.g., *"The runner now prefers..."*, *"The script automatically excludes..."*, *"The CI workflow now sets..."*).
- **Detail:** Mention specific flags, variables, timeouts, packages, or architectural details (e.g., specific package names like `numpy`, `netCDF4`, timeout values, directory paths).
- **Rationale:** Mention *why* the implementation choice was made (e.g., *"...ensuring reproducibility"*, *"...for easier debugging"*).
- **Mandatory File Citation:** Every single bullet point MUST end with the relative path of the file modified enclosed in parentheses, formatted as:
  `... description of the exact change and its rationale. (path/to/modified/file.ext)`
- Do not indent sub-bullets; keep all bullet points flat under their respective section header.

---

### 4. Guardrails & Prohibited Content
- **NO file lists or tables:** Do NOT include a "Files Changed", "Changed Files", "Walkthrough", "References", "Footnotes", or "Key Files" section or table at the bottom. The only place file paths should appear is inside the parentheses at the end of each bullet point.
- **NO reference links or citations:** Do not append link reference definitions (e.g., `[1]: ...` or `[file]: path`), issue links, or URLs at the end of the text.
- **NO generic headers:** Do not use markdown headers such as `### Changes`, `### Summary`, `## Overview`, or `### Details`. Use only the bold category headers specified in Section 2.
- **NO meta-announcements:** Output only the final pull request text. Do not include greetings, introductions (e.g., "Here is a draft of the pull request description:"), or closing conversational remarks.
