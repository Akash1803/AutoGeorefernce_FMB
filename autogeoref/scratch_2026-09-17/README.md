Working scripts and evidence from the 2026-09-17/18 sessions, copied out of the session scratchpad.
- autogeoref_demo.py: chain matcher (common_chain, rigid_fit, place with label-side, corroboration, rail and overlap rules), neighbour override merge.
- georef_village.py: village driver (seeds from hand files as poses, chain passes, block adjustment, status CSV, AUTOv2 files).
- block_adjust.py: whole-block rigid adjustment used for Thailavaram; manual_poses.py recovers rigid poses from a hand GeoJSON.
- neighbour_labels.py: glyph-based neighbour number reader (fragments; superseded by the two-reader transcription nb_override_35_04_077.json).
- make_nb_override.py: converts the read-fmb-neighbours workflow output (merged_077.json) into nb_override_<village>.json.
- topo_autodemo.py / load_demo.py / render_demo2.py / run_topo.py: QGIS-side review group, topology fix (railway conflicts unclipped) and rendering through the qgis MCP socket.
- spike_gcp/: Google tile export to GeoTIFF (google_sat_z20.xml, kizhi_sat_z20.tif), edge detection picture (kizhi_edges.png).
Village 35_04_077 transcription also lives in D:\Projects\Tambaram_Chengalpattu\FMB_Vector\35_04_077\neighbour_transcription\.
