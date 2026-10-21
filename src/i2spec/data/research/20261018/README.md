# Frozen research release: October 18, 2026

Four exact-J research models are installed here for explicit selection through
`i2spec.ResearchModel` and `i2spec research`. The production default remains
`i2spec2026o`. See [the API and scientific scope](../../../../../docs/research-models.md).

`manifest.json` records canonical IDs, role aliases, separate fitted/inherited/
calibration inventories, domain, exploratory benchmark scores, and SHA256 hashes
of every frozen source artifact. Files in each role directory are byte-for-byte
copies of the verified finalization-v2 snapshot; source paths in provenance are
historical records rather than runtime dependencies. The evaluator is the
maintained package module `i2spec.research`.

The four molecular inventories are 56+18, 61+11, 61+18 and 68+11, with four
calibration coefficients additional. No line/level offsets are fitted. Domain:
127I2, X v<=17, B v<=43, lower J<=200. Exact P(13) 43-0 exclusion retained.

`model.json` includes the staging-time `snapshot_prediction_replay_completed`
flag for exact provenance; the later successful `verification.json` and release
manifest report completed replay. Every integrated model reproduces the
60 absolute and 592 intraline reference predictions within 1 Hz. These are
reproduction checks, not independent predictive validation or uncertainty bounds.
