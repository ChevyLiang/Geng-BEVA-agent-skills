# Breeding preference rules

The functional-gene workflow has default breeding-direction rules and normally does not require a preference CSV.

Default interpretation:

- `Positive` or heading `promote`: Alternative/non-Nipponbare allele is favorable and Nipponbare is unfavorable.
- `Negative` or heading `delay`: Nipponbare allele is favorable and Alternative/non-Nipponbare is unfavorable.
- Heterozygous calls are reported separately from homozygous favorable/unfavorable calls.
- Earlier heading is the default breeding objective unless the user explicitly requests another objective.

When the user wants a different breeding objective, state that preference directly in the prompt. Explicit user preferences take precedence over the defaults for that analysis.
