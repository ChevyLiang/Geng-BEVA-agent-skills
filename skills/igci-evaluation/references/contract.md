# Common Python and manifest contract

Public analysis modules expose `run(input_path, out_dir="outputs", **kwargs) -> dict`.
Original named functions remain available for existing callers. Library functions raise on failure;
CLI processes then exit nonzero. The decorator writes an error manifest before re-raising.
Low-level helpers are not analysis modules and do not need a `run` function.

Required manifest fields: `schema_version`, `status`, `module`, `input`, `out_dir`, `outputs`, `warnings`, `metrics`.
`status` is `success` or `error`. `outputs` maps names to absolute paths of files actually written;
omit skipped artifacts. `warnings` is a list of strings; `metrics` is a JSON object.
Additional provenance may be placed in `extra` or preserved legacy fields.
`input` is the original requested file; preprocessors may use a cleaned derivative internally.

Each successful call writes `<module>.manifest.json` and `manifest.json` into `out_dir`.
`manifest.json` describes the latest call; use distinct run directories to retain independent results.
Invalid input is an error, not an empty successful analysis. An empty valid result is allowed.
No-informative-data numeric CSV entries are missing; JSON must not contain NaN/Infinity.
UTF-8 source and Markdown; CSV is UTF-8 with BOM for Excel compatibility.

Optional table/plot paths may also be present in existing `extra` provenance dictionaries.
Public outputs must not declare nonexistent files.
