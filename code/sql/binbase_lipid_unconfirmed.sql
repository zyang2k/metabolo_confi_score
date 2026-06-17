-- binbase_lipid_unconfirmed.sql — LIPIDOMICS (C18 reverse-phase) UNCONFIRMED de-orphan.
-- PURE SQL. Method below is C18 NEGATIVE; change 'negative' -> 'positive' for the pos run
-- (and the export filenames lipid_c18neg -> lipid_c18pos). Then locally:
--   python code/analysis/binbase_orphan_denoise_run.py --unconfirmed --tag lipid_c18neg

-- (1) distribution — confirm UNCONFIRMED count + that CONFIRMED has msms (paste back)
SELECT target_type, count(*) AS n,
       count(*) FILTER (WHERE msms IS NOT NULL AND msms<>'') AS with_msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | negative'
GROUP BY target_type ORDER BY n DESC;

-- (2) CONFIRMED bins (reference) -> EXPORT data/lipid_c18neg_bins.csv
SELECT id AS wiki_id, splash, version, accurate_mass AS precursor_mz,
       retention_time AS rt_sec, retention_index AS ri,
       pre_cursors_intensity AS precursor_intensity, name, adduct, ion_mode,
       fragment_of, fragmentation_parent_of, msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | negative'
  AND target_type = 'CONFIRMED' AND msms IS NOT NULL AND msms <> '';

-- (3) UNCONFIRMED candidate bins -> EXPORT data/lipid_c18neg_unconfirmed.csv
--     360k is heavy (export + run), so take a random 50k sample for a fast, statistically
--     tight fraction estimate. Drop "ORDER BY random() LIMIT 50000" for the full run.
SELECT id AS wiki_id, sample, splash, version, accurate_mass AS precursor_mz,
       retention_time AS rt_sec, retention_index AS ri,
       pre_cursors_intensity AS precursor_intensity, name, adduct, ion_mode,
       fragment_of, fragmentation_parent_of, msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | negative'
  AND target_type = 'UNCONFIRMED' AND msms IS NOT NULL AND msms <> ''
ORDER BY random() LIMIT 50000;
