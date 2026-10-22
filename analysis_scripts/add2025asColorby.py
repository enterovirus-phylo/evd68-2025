### Script to recolor the JSON to add 'study samples' as a colorby, so we can recolour them after adjusting the view for the manuscript figures

# After running Nextstrain tree, to create figures with the study samples
# more visible and coloured in, add a colorby specifically for the
# samples added on the date of the last run (update this below).
# This can viewed again in auspice, set to the correct zoom and size,
# and the phylogeny downloaded as SVG. The SVG can then be run through
# recolor_study_tips_manuFig.py to add the red outline to the study tips.

# run from root of repo (enterovirus/ev_githubs/enterovirus_d68-oct2025)

import json

IN_PATH  = "vp1/auspice/enterovirus_d68_vp1.json"
OUT_PATH = "vp1/auspice/enterovirus_d68_vp1_with_study_trait.json"
STUDY_DATE = "2026-10-10"

with open(IN_PATH) as f:
    data = json.load(f)

def annotate(node):
    if 'children' not in node:
        na = node.setdefault('node_attrs', {})
        is_study = na.get('date_added', {}).get('value') == STUDY_DATE
        na['study_sample'] = {"value": "This study (2025)" if is_study else "Background"}
    else:
        for c in node['children']:
            annotate(c)

annotate(data['tree'])

# Register as a coloring + filter so it's usable in Auspice's sidebar
colorings = data['meta'].setdefault('colorings', [])
if not any(c['key'] == 'study_sample' for c in colorings):
    colorings.append({
        "key": "study_sample",
        "title": "Study sequence",
        "type": "categorical",
        "scale": [
            ["This study (2025)", "#CC0000"],
            ["Background", "#AAAAAA"]
        ]
    })

filters = data['meta'].setdefault('filters', [])
if 'study_sample' not in filters:
    filters.append('study_sample')

with open(OUT_PATH, 'w') as f:
    json.dump(data, f)

def count_study(node, counts):
    if 'children' not in node:
        if node['node_attrs']['study_sample']['value'] == "This study (2025)":
            counts[0] += 1
    else:
        for c in node['children']:
            count_study(c, counts)

counts = [0]
count_study(data['tree'], counts)
print(f"Done. {counts[0]} study tips flagged out of full tree.")
