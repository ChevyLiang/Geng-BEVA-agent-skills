# PLINK runtime

This skill targets **genomic PLINK 1.9** for PCA, IBS, and IBD. PLINK 2 and PuTTY's `plink` are not substitutes for the current command wrappers.

## Windows-first deployment

For the primary Windows Codex workflow, put `plink.exe` at:

```text
skills/population-structure-analysis/tools/plink/windows-x64/plink.exe
```

The repository does not include the PLINK binary itself. The skill does not automatically install PLINK and does not write to `C:\Users\<user>\.geng-beva` or another user-profile cache.

Executable resolution order:

1. explicit `--plink-cmd`;
2. `PLINK_BIN` environment variable;
3. the skill-local executable above;
4. `plink` or `plink1.9` on `PATH`.

Each candidate is validated with `plink --version` and must identify itself as genomic PLINK 1.9. If none is available, the workflow reports the expected local path instead of attempting an installation.

For PLINK genome output, `DST` is used as the IBS sharing coefficient and `PI_HAT` as the estimated IBD proportion. IBD is not `1 - IBS`.
