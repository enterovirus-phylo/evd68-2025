"""
Recolor the outline (stroke) of 2025 study sample tips in an Auspice-exported SVG,
leaving fill (and everything else) untouched. To make study sample tips very visible.

Identity comes from the JSON's `study_sample` trait (added by add2025asColorby.py),
matched to SVG circle ids via Auspice's own id-sanitization rule:
any character outside [A-Za-z0-9_-] becomes '-'.

Example of use, running from base repo:

python3 scripts/recolor_study_tips_manuFig.py \
  --json vp1/auspice/enterovirus_d68_vp1_with_study_trait.json \
  --svg figs/b3.1_2025-08-17.svg \
  --out figs/recolored-b3.1_2025-08-17.svg

  python3 scripts/recolor_study_tips_manuFig.py \
  --json vp1/auspice/enterovirus_d68_vp1_with_study_trait.json \
  --svg figs/a2.3_2025-08-17.svg \
  --out figs/recolored-a2.3_2025-08-17.svg


"""
import json
import re
import argparse

def sanitize(name):
    return re.sub(r'[^A-Za-z0-9_-]', '-', name)

def get_study_names(json_path, flag_value="This study (2025)"):
    with open(json_path) as f:
        data = json.load(f)
    names = set()
    def walk(node):
        if 'children' not in node:
            val = node.get('node_attrs', {}).get('study_sample', {}).get('value')
            if val == flag_value:
                names.add(node['name'])
        else:
            for c in node['children']:
                walk(c)
    walk(data['tree'])
    return names

def recolor_svg(svg_path, study_names, out_path,
                 stroke_color="rgb(30, 30, 30)", #"rgb(204, 0, 0)", 
                 stroke_width="2.5"):
    with open(svg_path) as f:
        svg = f.read()

    sanitized_names = {sanitize(n): n for n in study_names}

    pattern = re.compile(r'<circle class="tip" id="tip_([^"]*)"([^>]*?)style="([^"]*)"([^/]*)/>')

    matched_ids = []
    unmatched_expected = set(sanitized_names.keys())

    def replace(m):
        raw_id, pre_attrs, style, post_attrs = m.groups()
        # match exact, or id with an appended "__accession" suffix
        is_study = raw_id in sanitized_names
        if not is_study:
            for san in sanitized_names:
                if raw_id.startswith(san + "__"):
                    is_study = True
                    unmatched_expected.discard(san)
                    break
        else:
            unmatched_expected.discard(raw_id)

        if not is_study:
            return m.group(0)

        matched_ids.append(raw_id)
        new_style = re.sub(r'stroke:\s*rgb\([^)]*\)', f'stroke: {stroke_color}', style)
        new_style = re.sub(r'stroke-width:\s*[0-9.]+', f'stroke-width: {stroke_width}', new_style)
        return f'<circle class="tip" id="tip_{raw_id}"{pre_attrs}style="{new_style}"{post_attrs}/>'

    new_svg = pattern.sub(replace, svg)

    with open(out_path, 'w') as f:
        f.write(new_svg)

    return matched_ids, unmatched_expected

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, help="Auspice JSON with study_sample trait")
    ap.add_argument("--svg", required=True, help="Exported SVG from Auspice")
    ap.add_argument("--out", required=True, help="Output SVG path")
    ap.add_argument("--flag-value", default="This study (2025)",
                     help="Value of study_sample marking a study tip")
    args = ap.parse_args()

    names = get_study_names(args.json, args.flag_value)
    matched, missing = recolor_svg(args.svg, names, args.out)
    print(f"Study samples in JSON: {len(names)}")
    print(f"Recolored in this SVG: {len(matched)}")
    if missing:
        print(f"Study samples not found in this SVG (likely just outside this zoom): {len(missing)}")