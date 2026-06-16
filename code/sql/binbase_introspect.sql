-- binbase_introspect.sql
-- Read-only introspection of carrot-prod `compound` to resolve the 4 unknowns
-- needed to run the orphan-denoise (reduce noise-ion UNCONFIRMED bins):
--   (1) exact column names for RT/RI, accurate mass, target-type, version, internalId
--   (2) msms serialization format
--   (3) RT/RI unit
--   (4) a HILIC method string + its CONFIRMED vs UNCONFIRMED counts
-- Nothing is written. Run all three and paste the output back.

-- ── A) full column inventory of `compound` (names + types) ──────────────────
SELECT column_name, data_type
FROM information_schema.columns
WHERE table_name = 'compound'
ORDER BY ordinal_position;

-- ── B) methods + row counts (pick a HILIC method from this list) ────────────
SELECT method, count(*) AS n
FROM compound
GROUP BY method
ORDER BY n DESC
LIMIT 50;

-- ── C) 3 sample rows so I can see msms format, RT/RI unit, and the
--        CONFIRMED/UNCONFIRMED target-type values in situ ────────────────────
--    (msms strings may be long — paste as-is, I'll parse them.)
SELECT *
FROM compound
LIMIT 3;
