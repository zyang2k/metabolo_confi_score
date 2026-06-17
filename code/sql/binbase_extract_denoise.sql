-- binbase_extract_denoise.sql  — PURE SQL (runs in any client; no psql meta-commands)
-- Read-only. Method: 5m hilic premier | orbitrap | beh amide | negative (a regular HILIC study).
--
-- HOW TO USE
--   Run each statement in your SQL tool. Paste back (1)'s small result.
--   EXPORT the result of (3) and (4) to CSV WITH A HEADER ROW, saved as:
--       (3) -> data/lcb_hilicneg_bins.csv
--       (4) -> data/lcb_hilicneg_orphans.csv
--   Then locally:  python code/analysis/binbase_orphan_denoise_run.py --from-csv
--   (the denoise = peak-level containment, can't be done in SQL, so it runs on the CSVs).

-- (1) target_type distribution — confirm CONFIRMED has msms; paste this back ------------
SELECT target_type,
       count(*)                                              AS n,
       count(*) FILTER (WHERE msms IS NOT NULL AND msms<>'') AS with_msms,
       count(*) FILTER (WHERE name LIKE 'unknown\_%')        AS unknown_named
FROM compound
WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
GROUP BY target_type
ORDER BY n DESC;

-- (2) top samples by floating-feature count (the largest is used as "one study") --------
SELECT sample, count(*) AS n_orphans
FROM compound
WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
  AND target_type <> 'CONFIRMED' AND msms IS NOT NULL AND msms <> ''
GROUP BY sample
ORDER BY n_orphans DESC
LIMIT 10;

-- (3) CONFIRMED bins (reference)  ->  EXPORT to data/lcb_hilicneg_bins.csv ---------------
--     Now carries CARROT's own ISF links (fragment_of / fragmentation_parent_of) for
--     validation #1 (labeled concordance).
SELECT id AS wiki_id, splash, version,
       accurate_mass    AS precursor_mz,
       retention_time   AS rt_sec,
       retention_index  AS ri,
       name, adduct, ion_mode,
       fragment_of, fragmentation_parent_of,
       msms
FROM compound
WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
  AND target_type = 'CONFIRMED' AND msms IS NOT NULL AND msms <> '';

-- (4) orphans across the top-50 injections  ->  EXPORT to data/lcb_hilicneg_orphans.csv
--     (bins from (3) are method-wide and reused as-is; only re-export this one.)
--     Carries fragment_of / fragmentation_parent_of = CARROT's OWN ISF links, for
--     concordance validation. Excludes ISTD / DELETED_BY_USER.
WITH samp AS (
  SELECT sample
  FROM compound
  WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
    AND target_type IN ('INVALID_TARGET','UNCONFIRMED') AND msms IS NOT NULL AND msms <> ''
  GROUP BY sample ORDER BY count(*) DESC LIMIT 50)
SELECT c.id AS wiki_id, c.sample, c.target_type, c.reason, c.splash,
       c.accurate_mass   AS precursor_mz,
       c.retention_time  AS rt_sec,
       c.retention_index AS ri,
       c.name, c.adduct, c.ion_mode,
       c.fragment_of, c.fragmentation_parent_of,
       c.msms
FROM compound c
WHERE c.method = '5m hilic premier | orbitrap | beh amide | negative'
  AND c.target_type IN ('INVALID_TARGET','UNCONFIRMED') AND c.msms IS NOT NULL AND c.msms <> ''
  AND c.sample IN (SELECT sample FROM samp);

-- ============================================================================
-- VALIDATION #2 — same-injection parent (sample_annotation_data)
-- ============================================================================
-- (5) introspect sample_annotation_data columns — paste back so I can fill the
--     sample-key and compound-key names into (6).
SELECT column_name, data_type
FROM information_schema.columns
WHERE table_name = 'sample_annotation_data'
ORDER BY ordinal_position;

-- (6) per-injection GENUINELY-DETECTED compounds for our 50 samples.
--     EXPORT to data/lcb_sample_detections.csv  (columns: sample, compound_id).
--     deleted=false AND replaced=false -> real detections, not gap-filled.
--     sample_annotation_data is a partitioned view: run EXPLAIN first; the
--     sample IN (...) filter should prune. If it scans everything, add an
--     `AND acquired >= 'YYYY-MM-01'` bound around this study's acquisition date.
WITH samp AS (
  SELECT sample FROM compound
  WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
    AND target_type IN ('INVALID_TARGET','UNCONFIRMED') AND msms IS NOT NULL AND msms <> ''
  GROUP BY sample ORDER BY count(*) DESC LIMIT 50)
SELECT sample, compound_id
FROM sample_annotation_data
WHERE sample IN (SELECT sample FROM samp)
  AND compound_id IS NOT NULL
  AND COALESCE(deleted, false)  = false
  AND COALESCE(replaced, false) = false;

-- ============================================================================
-- UNCONFIRMED candidate-bin audit (the "should this be promoted?" pile)
-- ============================================================================
-- (7) all UNCONFIRMED candidate bins for the method -> EXPORT to
--     data/lcb_hilicneg_unconfirmed.csv. Reference = confirmed bins from (3).
SELECT id AS wiki_id, sample, splash, version,
       accurate_mass   AS precursor_mz,
       retention_time  AS rt_sec,
       retention_index AS ri,
       name, adduct, ion_mode, fragment_of, fragmentation_parent_of, msms
FROM compound
WHERE method = '5m hilic premier | orbitrap | beh amide | negative'
  AND target_type = 'UNCONFIRMED' AND msms IS NOT NULL AND msms <> '';
