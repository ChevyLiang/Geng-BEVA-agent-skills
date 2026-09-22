# PLINK and external-process policy

The default functional-gene diagnostic workflow **does not require PLINK**. It streams the supplied VCF in Python and retains only the bundled RiceNavi-derived functional loci. This is intentional: it avoids unnecessary subprocess creation and is more robust on Windows, Linux, macOS, Codex sandboxes, and non-ASCII workspace paths.

`--plink-cmd` remains accepted for backward compatibility but is not needed for the normal workflow. Other Geng-BEVA skills (for example PCA/IBS/IBD) may still use PLINK, but PLINK is never auto-installed; population-structure analysis resolves the fixed skill-local executable or an explicitly configured existing installation.

## Windows/Codex process-launch fallback

If Codex itself cannot launch Python from the current workspace and reports `CreateProcessWithLogonW failed: 267`, treat this as a process working-directory problem rather than a VCF-format error. Do not rewrite the analysis into a temporary heredoc or piped `python -c` program. Instead:

1. resolve the Python interpreter, skill entry point, input VCF, and output directory to absolute paths;
2. use an ASCII-only working directory such as `%TEMP%\geng-beva-run`;
3. invoke the existing skill `main.py` from that working directory while keeping input/output paths absolute;
4. verify declared outputs and the requested sample ID after completion.
