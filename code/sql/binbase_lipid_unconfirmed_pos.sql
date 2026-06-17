-- binbase_lipid_unconfirmed_pos.sql — LIPIDOMICS C18 POSITIVE UNCONFIRMED de-orphan.
-- PURE SQL. Then locally: python code/analysis/binbase_orphan_denoise_run.py --unconfirmed --tag lipid_c18pos

-- (1) distribution (paste back)
SELECT target_type, count(*) AS n,
       count(*) FILTER (WHERE msms IS NOT NULL AND msms<>'') AS with_msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | positive'
GROUP BY target_type ORDER BY n DESC;

-- (2) CONFIRMED bins (reference) -> EXPORT data/lipid_c18pos_bins.csv
SELECT id AS wiki_id, splash, version, accurate_mass AS precursor_mz,
       retention_time AS rt_sec, retention_index AS ri,
       pre_cursors_intensity AS precursor_intensity, name, adduct, ion_mode,
       fragment_of, fragmentation_parent_of, msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | positive'
  AND target_type = 'CONFIRMED' AND msms IS NOT NULL AND msms <> '';

-- (3) UNCONFIRMED candidates (random 50k sample for speed) -> EXPORT data/lipid_c18pos_unconfirmed.csv
--     drop "ORDER BY random() LIMIT 50000" for the full run.
SELECT id AS wiki_id, sample, splash, version, accurate_mass AS precursor_mz,
       retention_time AS rt_sec, retention_index AS ri,
       pre_cursors_intensity AS precursor_intensity, name, adduct, ion_mode,
       fragment_of, fragmentation_parent_of, msms
FROM compound
WHERE method = '5m splash one premier | orbitrap | beh c18 | positive'
  AND target_type = 'UNCONFIRMED' AND msms IS NOT NULL AND msms <> ''
ORDER BY random() LIMIT 50000;
